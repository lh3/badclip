//! The `merge` subcommand: collapse per-read breakend records (the output of
//! `extract`) into consensus structural-variant calls with read counts.
//!
//! This is a simplified reimplementation of `gc_cmd_merge` in `test/minisv.js`.
//! badclip only has the breakend type — no INS/DEL/DUP/INV typing, no
//! `SVTYPE`/`SVLEN`, no per-sample `source`, and no centromere or RT (TSD/polyA)
//! annotation — so minisv's SVTYPE/SVLEN clustering constraints, the `-d` allele
//! length filter, the centromere filter (`-e`), and the whole RT branch
//! (`-r`/`-R`) are dropped.
//!
//! Unlike minisv (which assumes an upstream `sort -k1,1 -k2,2n`), this loads all
//! records into memory and sorts them itself; the extract output is not huge.
//!
//! Algorithm: sort records by `(ctg, pos)`, then sweep a window of active
//! clusters. Each incoming record joins the active cluster whose members it most
//! often matches (`same_sv`), else starts a new cluster. When a cluster falls
//! out of the window it is emitted if it passes the read-count filters.

use std::collections::{BTreeMap, VecDeque};
use std::io::{self, BufRead, Write};

use crate::io::open_reader;

/// Options for `merge`, mirroring minisv's `-c -s -w -A` (plus `-M` for the
/// per-cluster compare cap and the badclip-specific `-m`/`-C`/`-S` filters).
pub struct MergeOpts {
    /// `extract` output (gzip ok; "-" for stdin).
    pub input: String,
    /// Minimum read count to emit a call.
    pub min_cnt: i64,
    /// Minimum read count on each strand.
    pub min_cnt_strand: i64,
    /// Clustering window size (bp).
    pub win_size: i64,
    /// Cap on active clusters (flush trigger).
    pub max_allele: i64,
    /// Maximum members compared per cluster (deterministic cap).
    pub max_check: i64,
    /// Drop input breakends whose eseq quality (`equal` tag) is below this.
    pub min_equal: i64,
    /// Drop input breakends whose col-7 mapq (join `max`, clip mapq) is below this.
    pub min_mapq: i64,
    /// Drop input **join** breakends whose smaller per-hit mapq is below this
    /// (no effect on clips).
    pub min_mapq_min: i64,
    /// Treat input as combined `merge` output (sample-merge) rather than
    /// `extract` output. Changes only input parsing and output INFO emission.
    pub merge_input: bool,
    /// (-m only) drop input lines whose total count (Σ fwd+rev) is below this.
    pub min_cnt_in: i64,
    /// (-m only) drop input lines whose per-strand count is below this.
    pub min_cnt_strand_in: i64,
}

/// One parsed `extract` breakend record. Every extract record is a breakend, so
/// there is no indel variant. Only fields needed for clustering/emission are
/// kept (the col-8 INFO is not — merge INFO is derived from scratch).
struct Rec {
    ctg: String,       // col 0
    pos: i64,          // col 1 (0-based)
    ori: [u8; 2],      // col 2, two chars each b'>'/b'<'/b'.'
    ctg2: String,      // col 3 ("." for a clip)
    pos2: Option<i64>, // col 4 (None when ".")
    // Per-hit mapqs from the `mapq=` INFO tag, already aligned to the output
    // endpoints (mapq1 -> ctg1:pos1, mapq2 -> ctg2:pos2; extract swaps them on
    // flip). `mapq2` is 0 for a clip. Fall back to (col-7 mapq, 0) if the tag is
    // absent. Used for the per-end `avg_mapq=q1,q2`.
    mapq1: i64,
    mapq2: i64,
    strand: u8, // col 7, b'+' or b'-'
    name: String,      // col 5, qname
    source: String,    // `source=` INFO tag (col 8); "." if absent
    // Read counts carried by a `-m` (sample-merge) input line: `fwd`/`rev` summed
    // over the `count=` tag's sources, `fr`/`rf` from `count_fr=`/`count_rf=`. All
    // 0 in extract mode, where write_sv counts members by strand/ori instead.
    fwd: i64,
    rev: i64,
    fr: i64,
    rf: i64,
    // Source labels parsed from the `count=` tag of a `-m` input line (the `src:`
    // prefixes). Empty in extract mode and for label-less (bare `f,r`) counts.
    // Unioned across a cluster to give the col-7 source count.
    sources: Vec<String>,
}

/// An open cluster of records deemed to describe the same breakend.
struct Cluster {
    ctg: String,
    pos_max: i64,
    members: Vec<Rec>,
}

/// Parse one `extract` line into a `Rec`; `None` for malformed lines (`< 9`
/// fields), lines whose col-7 mapq is below `min_mapq`, lines whose `equal`
/// quality is below `min_equal` (lines lacking an `equal` tag are kept), or
/// **join** lines whose smaller per-hit mapq (from the `mapq=` INFO tag) is below
/// `min_mapq_min` (clips are exempt), consistent with the other subcommands.
fn parse_rec(line: &str, min_equal: i64, min_mapq: i64, min_mapq_min: i64) -> Option<Rec> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 9 {
        return None;
    }
    let ob = f[2].as_bytes();
    if ob.len() < 2 {
        return None;
    }
    // Drop breakends whose col-7 mapq (the join's max, or a clip's mapq) is below
    // the threshold.
    let mapq: i64 = f[6].parse().ok()?;
    if mapq < min_mapq {
        return None;
    }
    // `.` (a clip's missing mate) fails to parse -> None, the sentinel `same_sv`
    // wants; also gates the -p clip exemption below.
    let pos2 = f[4].parse::<i64>().ok();
    // Per-hit mapq pair from the `mapq=` INFO tag (already aligned to the output
    // endpoints). Fall back to (col-7 mapq, 0) when the tag is absent — col-7 is
    // the clip's own mapq, so this matches a clip's (mapq1, 0).
    let (mapq1, mapq2) = f[8]
        .split(';')
        .find_map(|kv| kv.strip_prefix("mapq="))
        .and_then(|v| {
            let mut it = v.split(',');
            Some((it.next()?.parse().ok()?, it.next()?.parse().ok()?))
        })
        .unwrap_or((mapq, 0));
    // -p: drop join breakends whose *smaller* per-hit mapq is below the threshold.
    // Clips (pos2 absent) are exempt: a clip's mapq2 is 0, so its min is always 0.
    if pos2.is_some() && mapq1.min(mapq2) < min_mapq_min {
        return None;
    }
    // Drop low-quality-eseq breakends; keep those without an `equal` tag.
    if let Some(q) = f[8]
        .split(';')
        .find_map(|kv| kv.strip_prefix("equal="))
        .and_then(|v| v.parse::<i64>().ok())
        && q < min_equal
    {
        return None;
    }
    Some(Rec {
        ctg: f[0].to_string(),
        pos: f[1].parse().ok()?,
        ori: [ob[0], ob[1]],
        ctg2: f[3].to_string(),
        pos2,
        name: f[5].to_string(),
        mapq1,
        mapq2,
        strand: f[7].bytes().next()?,
        // Dataset label from the `source=` INFO tag (order-independent); "." when
        // absent, though extract always emits it.
        source: f[8]
            .split(';')
            .find_map(|kv| kv.strip_prefix("source="))
            .unwrap_or(".")
            .to_string(),
        // Not used in extract mode (write_sv counts members by strand/ori).
        fwd: 0,
        rev: 0,
        fr: 0,
        rf: 0,
        sources: Vec::new(),
    })
}

/// Parse a `count=` INFO value: sum its `|`-separated `[src:]f,r` entries into
/// `(fwd, rev)` and collect the `src:` labels. The source-prefix is optional (an
/// entry may be a bare `f,r`, as in `-m` output) — such entries contribute to the
/// counts but no label. Mirrors minisv's `gc_get_count_gsv`.
fn parse_count(info: &str) -> (i64, i64, Vec<String>) {
    let mut fwd = 0;
    let mut rev = 0;
    let mut sources = Vec::new();
    if let Some(v) = info.split(';').find_map(|kv| kv.strip_prefix("count=")) {
        for part in v.split('|') {
            // Split off an optional "src:" prefix (source names carry no ':').
            let (label, fr) = match part.rsplit_once(':') {
                Some((lbl, fr)) => (Some(lbl), fr),
                None => (None, part),
            };
            let mut it = fr.split(',');
            if let (Some(a), Some(b)) = (it.next(), it.next())
                && let (Ok(a), Ok(b)) = (a.parse::<i64>(), b.parse::<i64>())
            {
                fwd += a;
                rev += b;
                if let Some(lbl) = label {
                    sources.push(lbl.to_string());
                }
            }
        }
    }
    (fwd, rev, sources)
}

/// Parse one integer-valued INFO tag (`key=N`); 0 if absent/unparseable.
fn info_int(info: &str, key: &str) -> i64 {
    info.split(';')
        .find_map(|kv| kv.strip_prefix(key))
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

/// Parse one combined-`merge`-output line into a `Rec` for `-m` (sample-merge)
/// mode. Same 9 columns as `extract`, but the INFO fields differ: mapq comes from
/// `avg_mapq=q1,q2` (col f[6] is a count here, not mapq) and the per-strand read
/// counts from `count=` (summed over sources). `None` for malformed lines, lines
/// failing the per-line `-C`/`-S` count filters, or lines failing `-q`/`-p` on the
/// avg-mapq-derived q1/q2 (`-Q`/`equal` does not apply — merge output has none).
fn parse_rec_merge(line: &str, opts: &MergeOpts) -> Option<Rec> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 9 {
        return None;
    }
    let ob = f[2].as_bytes();
    if ob.len() < 2 {
        return None;
    }
    let pos2 = f[4].parse::<i64>().ok();
    // mapq from `avg_mapq=q1,q2` (0,0 if absent).
    let (mapq1, mapq2) = f[8]
        .split(';')
        .find_map(|kv| kv.strip_prefix("avg_mapq="))
        .and_then(|v| {
            let mut it = v.split(',');
            Some((it.next()?.parse().ok()?, it.next()?.parse().ok()?))
        })
        .unwrap_or((0, 0));
    // -q: drop when the larger per-end avg mapq is below the threshold.
    if mapq1.max(mapq2) < opts.min_mapq {
        return None;
    }
    // -p: drop joins whose smaller per-end avg mapq is below the threshold
    // (clips, pos2 absent, are exempt — their q2 is 0).
    if pos2.is_some() && mapq1.min(mapq2) < opts.min_mapq_min {
        return None;
    }
    let (fwd, rev, sources) = parse_count(f[8]);
    // Per-line -C (total) and -S (per-strand) filters.
    if fwd + rev < opts.min_cnt_in {
        return None;
    }
    if fwd < opts.min_cnt_strand_in || rev < opts.min_cnt_strand_in {
        return None;
    }
    Some(Rec {
        ctg: f[0].to_string(),
        pos: f[1].parse().ok()?,
        ori: [ob[0], ob[1]],
        ctg2: f[3].to_string(),
        pos2,
        name: ".".to_string(),
        mapq1,
        mapq2,
        strand: f[7].bytes().next()?,
        source: ".".to_string(),
        fwd,
        rev,
        fr: info_int(f[8], "count_fr="),
        rf: info_int(f[8], "count_rf="),
        sources,
    })
}

/// Whether two records describe the same breakend (simplified from minisv:
/// no SVTYPE/SVLEN checks).
fn same_sv(v: &Rec, w: &Rec, win: i64) -> bool {
    if v.ctg != w.ctg {
        return false;
    }
    // ori must match, except the inversion pair "><" <-> "<>" is compatible.
    if v.ori != w.ori {
        let inv = (v.ori == *b"><" && w.ori == *b"<>") || (v.ori == *b"<>" && w.ori == *b"><");
        if !inv {
            return false;
        }
    }
    if (v.pos - w.pos).abs() > win {
        return false;
    }
    if v.ctg2 != w.ctg2 {
        return false;
    }
    // pos2 window: minisv relies on `NaN > x == false`, so the numeric check is
    // skipped when either side is a clip (pos2 absent). Two clips then match on
    // ctg2 equality alone; a clip vs a join is already separated by ctg2.
    if let (Some(a), Some(b)) = (v.pos2, w.pos2)
        && (a - b).abs() > win
    {
        return false;
    }
    true
}

/// Apply the read-count filters and, if the cluster passes, emit one merged line.
fn write_sv(out: &mut impl Write, cl: &Cluster, opts: &MergeOpts) -> io::Result<()> {
    let s = &cl.members;
    if s.is_empty() {
        return Ok(());
    }
    let v = &s[s.len() / 2]; // representative (upper-middle; stable under sort)

    if opts.merge_input {
        return write_sv_merge(out, s, v, opts);
    }

    // Per-end mapq sums. mapq1 tracks ctg1:pos1, mapq2 tracks ctg2:pos2; both are
    // already output-endpoint-aligned per record (extract swaps the pair on flip),
    // and `same_sv` clusters only records with matching pos1 *and* pos2, so every
    // member's end1/end2 sit at the same loci — including the "><"/"<>" inversion
    // pair, whose coordinate-based canonicalization keeps end1 = the lower locus.
    let (mut mapq1_sum, mut mapq2_sum) = (0i64, 0i64);
    let mut cnt = [0i64, 0i64]; // [+, -] cluster totals (for the filter)
    let (mut fr, mut rf) = (0i64, 0i64); // ori "><" and "<>"
    // Per-source [+, -] counts and read-name lists; BTreeMap keeps sources
    // alphabetical in output. Read names stay in cluster-member order per source.
    let mut per_source: BTreeMap<&str, [i64; 2]> = BTreeMap::new();
    let mut reads_by_source: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for m in s {
        mapq1_sum += m.mapq1;
        mapq2_sum += m.mapq2;
        let strand_ix = if m.strand == b'+' { 0 } else { 1 };
        cnt[strand_ix] += 1;
        per_source.entry(&m.source).or_default()[strand_ix] += 1;
        reads_by_source.entry(&m.source).or_default().push(&m.name);
        if m.ori == *b"><" {
            fr += 1;
        } else if m.ori == *b"<>" {
            rf += 1;
        }
    }

    // Count filter (minisv's RT relaxation branch is dropped, so unconditional).
    if (s.len() as i64) < opts.min_cnt {
        return Ok(());
    }
    if cnt[0] < opts.min_cnt_strand || cnt[1] < opts.min_cnt_strand {
        return Ok(());
    }

    // f64::round is half-away-from-zero, matching JS toFixed(0). Per-end: q1 for
    // the ctg1:pos1 side, q2 for the ctg2:pos2 side (q2 is 0 for a clip cluster).
    let n = s.len() as f64;
    let q1 = (mapq1_sum as f64 / n).round() as i64;
    let q2 = (mapq2_sum as f64 / n).round() as i64;

    // count=src:f,r|... over sources present in the cluster, alphabetical.
    let count = per_source
        .iter()
        .map(|(src, [f, r])| format!("{src}:{f},{r}"))
        .collect::<Vec<_>>()
        .join("|");
    let mut info = format!("avg_mapq={q1},{q2};count={count}");
    if fr + rf > 0 {
        info.push_str(&format!(";count_fr={fr};count_rf={rf}"));
        if fr * rf == 0 && v.ctg == v.ctg2 {
            info.push_str(";foldback");
        }
    }
    // reads=src:name,...|... — same sources/order as count=, names per source.
    let reads = reads_by_source
        .iter()
        .map(|(src, ns)| format!("{src}:{}", ns.join(",")))
        .collect::<Vec<_>>()
        .join("|");
    info.push_str(&format!(";reads={reads}"));

    let pos2 = v.pos2.map_or_else(|| ".".to_string(), |p| p.to_string());
    let ori = std::str::from_utf8(&v.ori).unwrap_or(".");
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}",
        v.ctg,
        v.pos,
        ori,
        v.ctg2,
        pos2,
        s.len(),
        v.strand as char,
        info,
    )
}

/// Emit one sample-merged line for `-m` mode. Aggregates the members' `count=`
/// read totals (per strand) and `count_fr`/`count_rf`, filters on the read-based
/// `-c`/`-s`, and derives INFO from scratch — `avg_mapq=q1,q2` (mean of members'
/// per-end avg mapq), `count=F,R` (no per-source breakdown), and, for inversions,
/// `count_fr`/`count_rf`/`foldback`. No `reads=` tag. The col-7 count slot is the
/// number of distinct **sources** (samples) in the cluster. Output is itself
/// `-m`-parseable, so merges can be chained (though it carries no source labels,
/// so a chained pass falls back to member count for col-7).
fn write_sv_merge(
    out: &mut impl Write,
    s: &[Rec],
    v: &Rec,
    opts: &MergeOpts,
) -> io::Result<()> {
    let (mut mapq1_sum, mut mapq2_sum) = (0i64, 0i64);
    let (mut fwd, mut rev, mut fr, mut rf) = (0i64, 0i64, 0i64, 0i64);
    let mut src_set: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for m in s {
        mapq1_sum += m.mapq1;
        mapq2_sum += m.mapq2;
        fwd += m.fwd;
        rev += m.rev;
        fr += m.fr;
        rf += m.rf;
        for src in &m.sources {
            src_set.insert(src.as_str());
        }
    }
    // col-7 = distinct source count; fall back to member count when the input
    // carried no source labels (e.g. a chained `-m` pass over `count=F,R`).
    let n_src = if src_set.is_empty() {
        s.len()
    } else {
        src_set.len()
    };

    // Read-based -c (total) and -s (per-strand) filters.
    if fwd + rev < opts.min_cnt {
        return Ok(());
    }
    if fwd < opts.min_cnt_strand || rev < opts.min_cnt_strand {
        return Ok(());
    }

    let n = s.len() as f64;
    let q1 = (mapq1_sum as f64 / n).round() as i64;
    let q2 = (mapq2_sum as f64 / n).round() as i64;
    let mut info = format!("avg_mapq={q1},{q2};count={fwd},{rev}");
    if fr + rf > 0 {
        info.push_str(&format!(";count_fr={fr};count_rf={rf}"));
        if fr * rf == 0 && v.ctg == v.ctg2 {
            info.push_str(";foldback");
        }
    }

    let pos2 = v.pos2.map_or_else(|| ".".to_string(), |p| p.to_string());
    let ori = std::str::from_utf8(&v.ori).unwrap_or(".");
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}",
        v.ctg,
        v.pos,
        ori,
        v.ctg2,
        pos2,
        n_src,
        v.strand as char,
        info,
    )
}

/// Read `extract` output, cluster/merge breakends, and print consensus calls.
pub fn run(opts: &MergeOpts) -> io::Result<()> {
    // Load everything, then sort by (ctg, pos). Rust `str` ordering equals
    // `LC_ALL=C sort` byte order (the pipeline convention); the sort is stable,
    // so equal-key records keep input order and the representative pick is
    // deterministic.
    let mut recs: Vec<Rec> = Vec::new();
    for line in open_reader(&opts.input)?.lines() {
        let line = line?;
        // `-m`: input is combined `merge` output; otherwise `extract` output.
        let rec = if opts.merge_input {
            parse_rec_merge(&line, opts)
        } else {
            parse_rec(&line, opts.min_equal, opts.min_mapq, opts.min_mapq_min)
        };
        if let Some(r) = rec {
            recs.push(r);
        }
    }
    recs.sort_by(|a, b| a.ctg.cmp(&b.ctg).then(a.pos.cmp(&b.pos)));

    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    let max_check = opts.max_check.max(0) as usize;
    let mut active: VecDeque<Cluster> = VecDeque::new();
    for v in recs {
        // Flush clusters that can no longer accept v. `active.len()` is checked
        // before pushing v, matching minisv's pre-increment semantics.
        while let Some(front) = active.front() {
            if front.ctg != v.ctg
                || v.pos - front.pos_max > opts.win_size
                || active.len() as i64 > opts.max_allele
            {
                write_sv(&mut out, &active.pop_front().unwrap(), opts)?;
            } else {
                break;
            }
        }

        // Assign v to the cluster whose members it most often matches. Only the
        // first `max_check` members are compared (a deterministic cap replacing
        // minisv's reservoir sampling), so counts on very deep clusters are
        // capped rather than scaled.
        let mut best_i: Option<usize> = None;
        let mut best_c = 0;
        for (i, cl) in active.iter().enumerate() {
            let cap = cl.members.len().min(max_check);
            let mut c = 0;
            for m in &cl.members[..cap] {
                if same_sv(m, &v, opts.win_size) {
                    c += 1;
                }
            }
            if c > best_c {
                best_c = c;
                best_i = Some(i);
            }
        }
        match best_i {
            Some(i) => {
                let cl = &mut active[i];
                cl.pos_max = cl.pos_max.max(v.pos);
                cl.members.push(v);
            }
            None => active.push_back(Cluster {
                ctg: v.ctg.clone(),
                pos_max: v.pos,
                members: vec![v],
            }),
        }
    }
    while let Some(cl) = active.pop_front() {
        write_sv(&mut out, &cl, opts)?;
    }
    Ok(())
}
