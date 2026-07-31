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

/// Options for `merge`, mirroring minisv's `-c -s -w -A -C`.
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
    mapq: i64,         // col 6
    strand: u8,        // col 7, b'+' or b'-'
    name: String,      // col 5, qname
    source: String,    // `source=` INFO tag (col 8); "." if absent
}

/// An open cluster of records deemed to describe the same breakend.
struct Cluster {
    ctg: String,
    pos_max: i64,
    members: Vec<Rec>,
}

/// Parse one `extract` line into a `Rec`; `None` for malformed lines (`< 9`
/// fields), consistent with the other subcommands.
fn parse_rec(line: &str) -> Option<Rec> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 9 {
        return None;
    }
    let ob = f[2].as_bytes();
    if ob.len() < 2 {
        return None;
    }
    Some(Rec {
        ctg: f[0].to_string(),
        pos: f[1].parse().ok()?,
        ori: [ob[0], ob[1]],
        ctg2: f[3].to_string(),
        // "." (a clip's missing mate) fails to parse -> None, which is exactly
        // the sentinel we want (see `same_sv`).
        pos2: f[4].parse().ok(),
        name: f[5].to_string(),
        mapq: f[6].parse().ok()?,
        strand: f[7].bytes().next()?,
        // Dataset label from the `source=` INFO tag (order-independent); "." when
        // absent, though extract always emits it.
        source: f[8]
            .split(';')
            .find_map(|kv| kv.strip_prefix("source="))
            .unwrap_or(".")
            .to_string(),
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

    let mut mapq_sum = 0i64;
    let mut cnt = [0i64, 0i64]; // [+, -] cluster totals (for the filter)
    let (mut fr, mut rf) = (0i64, 0i64); // ori "><" and "<>"
    let mut names: Vec<&str> = Vec::with_capacity(s.len());
    // Per-source [+, -] counts; BTreeMap keeps sources alphabetical in output.
    let mut per_source: BTreeMap<&str, [i64; 2]> = BTreeMap::new();
    for m in s {
        mapq_sum += m.mapq;
        let strand_ix = if m.strand == b'+' { 0 } else { 1 };
        cnt[strand_ix] += 1;
        per_source.entry(&m.source).or_default()[strand_ix] += 1;
        if m.ori == *b"><" {
            fr += 1;
        } else if m.ori == *b"<>" {
            rf += 1;
        }
        names.push(&m.name);
    }

    // Count filter (minisv's RT relaxation branch is dropped, so unconditional).
    if (s.len() as i64) < opts.min_cnt {
        return Ok(());
    }
    if cnt[0] < opts.min_cnt_strand || cnt[1] < opts.min_cnt_strand {
        return Ok(());
    }

    // f64::round is half-away-from-zero, matching JS toFixed(0).
    let avg_mapq = (mapq_sum as f64 / s.len() as f64).round() as i64;

    // count=src:f,r|... over sources present in the cluster, alphabetical.
    let count = per_source
        .iter()
        .map(|(src, [f, r])| format!("{src}:{f},{r}"))
        .collect::<Vec<_>>()
        .join("|");
    let mut info = format!("avg_mapq={avg_mapq};count={count}");
    if fr + rf > 0 {
        info.push_str(&format!(";count_fr={fr};count_rf={rf}"));
        if fr * rf == 0 && v.ctg == v.ctg2 {
            info.push_str(";foldback");
        }
    }
    info.push_str(&format!(";reads={}", names.join(",")));

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

/// Read `extract` output, cluster/merge breakends, and print consensus calls.
pub fn run(opts: &MergeOpts) -> io::Result<()> {
    // Load everything, then sort by (ctg, pos). Rust `str` ordering equals
    // `LC_ALL=C sort` byte order (the pipeline convention); the sort is stable,
    // so equal-key records keep input order and the representative pick is
    // deterministic.
    let mut recs: Vec<Rec> = Vec::new();
    for line in open_reader(&opts.input)?.lines() {
        if let Some(r) = parse_rec(&line?) {
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
