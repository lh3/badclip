//! The `extract` subcommand: extract breakends from read alignments.
//!
//! A read may map as several chimeric hits. Sorted along the read, each read
//! end that is soft-clipped away from an alignment produces a *trim* breakend
//! (connected to nothing; only with `-T`), each junction between two adjacent
//! hits produces a *split* breakend, and each long indel (`-g`) inside a hit's
//! CIGAR produces a *del* or *ins* breakend. All records share one 9-column,
//! TAB-delimited layout:
//!
//! ```text
//! ctg1  pos1  ori  ctg2  pos2  qname  mapq  strand  source=..;ctype=..;idx=..;n_aln=..;aln_len=..;qlen=..;mapq=..
//! ```
//!
//! where `ori` is two characters, each `>`/`<`, or `.` for a missing mate. The
//! INFO column opens with `source=<name>` (the `-s` dataset label; first only
//! for readability — tag order is not significant), `ctype` (`trim`, `split`,
//! `del` or `ins`), `idx`, and `n_aln` (the read's alignment count that passed
//! filtering), then three tags:
//! - `aln_len=len1,len2` — alignment block lengths (PAF col 11) of the hits at
//!   `ctg1:pos1` and `ctg2:pos2`; for a trim `len2` is `0`.
//! - `qlen=left,middle,right` — query lengths whose sum is the read length.
//!   `left` is on the `ctg1:pos1` side, `right` on the `ctg2:pos2` side. For a
//!   trim, `right` is the clipped length and `middle` is `0`; for a split,
//!   `middle` is the (possibly negative) query gap between the two hits.
//! - `mapq=mapq1,mapq2` — the mapq of each hit, `mapq1` for `ctg1:pos1` and
//!   `mapq2` for `ctg2:pos2`; for a trim `mapq2` is `0`. (Distinct from the
//!   field-7 mapq column, which is the `max` of a split's two mapqs.)
//!
//! A del/ins breakend is written like the split an aligner would report for it
//! (`ctg t >> ctg t+len` for a deletion), with both "hits" being the one
//! alignment. An insertion is written `ctg L >> ctg R`, its leftmost/rightmost
//! equivalent placements, and adds `inslen` plus a second eseq window
//! (`elen2`/`equal2`/`eseq2`) for its second breakpoint on the read.

use std::io::{self, BufRead, Write};

use crate::io::open_reader;
use crate::paf::{Gap, Hit, Strand, parse_paf_line, parse_paf_unmapped};

/// Options for `badclip extract`.
pub struct ExtractOpts {
    /// Input path; `"-"` means stdin.
    pub input: String,
    /// Read PAF instead of an alignment file (BAM/CRAM/SAM).
    pub paf: bool,
    /// Reference FASTA (faidx-indexed), required for CRAM input.
    pub reference: Option<String>,
    /// File listing ALT contig names (alignment-file input); `None` = none.
    pub alt: Option<String>,
    /// Dataset name, emitted as the `source=` INFO tag on every record.
    pub source: String,
    /// Minimum clip length to report a trim breakend.
    pub min_clip: i64,
    /// Emit trim breakends (read ends clipped away from every alignment).
    pub trim: bool,
    /// Minimum CIGAR indel length for a del/ins breakend (`<= 0` disables).
    pub min_gap: i64,
    /// Drop hits with mapq below this value (0 = keep everything).
    pub min_mapq: i64,
    /// Drop hits whose alignment block length is below this (0 = keep all).
    pub min_aln_len: i64,
    /// Flanking read sequence extracted on each side of a breakend (BAM only).
    pub flank: i64,
    /// Maximum extracted sequence length; longer windows omit `eseq` (BAM only).
    pub max_eseq: i64,
    /// Drop breakends whose eseq quality (`equal`) is below this (0 = keep all).
    pub min_equal: i64,
}

/// Read lengths up to this are tallied in a fixed counting array; only longer
/// reads are kept individually, so a billion-short-read BAM needs O(1) memory.
const MAX_COUNTED_LEN: i64 = 65535;

/// Per-read input statistics, printed to stderr at the end of a run.
///
/// A "read" is one primary record (alignment input; secondary/supplementary
/// records are ignored, unmapped ones are counted and also tallied in
/// `n_unmapped`) or one PAF read group; reads are counted before the `-a`
/// filter, so a fully filtered read still counts. `primary_bases` sums the query span (`qe - qs`) of each read's
/// primary alignment — PAF carries no primary flag (all kept lines are
/// `tp:A:P`), so the read's longest query span stands in for it there.
pub(crate) struct Stats {
    n_reads: u64,
    n_unmapped: u64,
    primary_bases: u64,
    total_len: u64,
    /// `short_cnt[l]` = number of reads of length `l <= MAX_COUNTED_LEN`.
    short_cnt: Vec<u64>,
    /// Individual lengths of reads longer than `MAX_COUNTED_LEN` (rare).
    long_lens: Vec<i64>,
}

impl Default for Stats {
    fn default() -> Self {
        Stats {
            n_reads: 0,
            n_unmapped: 0,
            primary_bases: 0,
            total_len: 0,
            short_cnt: vec![0; (MAX_COUNTED_LEN + 1) as usize],
            long_lens: Vec::new(),
        }
    }
}

impl Stats {
    /// Account one read: its length and its primary alignment's query span.
    pub(crate) fn add_read(&mut self, qlen: i64, primary_span: i64) {
        self.n_reads += 1;
        self.primary_bases += primary_span.max(0) as u64;
        let qlen = qlen.max(0);
        self.total_len += qlen as u64;
        if qlen <= MAX_COUNTED_LEN {
            self.short_cnt[qlen as usize] += 1;
        } else {
            self.long_lens.push(qlen);
        }
    }

    /// Account one unmapped read: it counts toward the read total and the N50
    /// but contributes no primary-aligned bases.
    pub(crate) fn add_unmapped(&mut self, qlen: i64) {
        self.n_unmapped += 1;
        self.add_read(qlen, 0);
    }

    /// N50 read length: walking reads from longest to shortest, the length at
    /// which the running sum first reaches half the total read length. Long
    /// reads are walked individually, then the counting array from the top.
    fn n50(&mut self) -> i64 {
        let total = self.total_len;
        let mut cum = 0u64;
        self.long_lens.sort_unstable_by(|a, b| b.cmp(a));
        for &len in &self.long_lens {
            cum += len as u64;
            if cum * 2 >= total {
                return len;
            }
        }
        for len in (1..=MAX_COUNTED_LEN).rev() {
            let cnt = self.short_cnt[len as usize];
            if cnt == 0 {
                continue;
            }
            cum += cnt * len as u64;
            if cum * 2 >= total {
                // The threshold falls inside this length bin, possibly short of
                // its last copy; every copy has the same length, so return it.
                return len;
            }
        }
        0
    }

    /// Print the summary to stderr.
    fn report(&mut self) {
        eprintln!("number of reads: {}", self.n_reads);
        eprintln!("number of unmapped reads: {}", self.n_unmapped);
        eprintln!("total bases in primary alignments: {}", self.primary_bases);
        eprintln!("read N50: {}", self.n50());
    }
}

#[cfg(test)]
mod tests {
    use super::{MAX_COUNTED_LEN, Stats};

    /// N50 computed the straightforward way, for cross-checking.
    fn naive_n50(lens: &[i64]) -> i64 {
        let mut lens = lens.to_vec();
        lens.sort_unstable_by(|a, b| b.cmp(a));
        let total: i64 = lens.iter().sum();
        let mut cum = 0;
        for &l in &lens {
            cum += l;
            if cum * 2 >= total {
                return l;
            }
        }
        0
    }

    fn stats_n50(lens: &[i64]) -> i64 {
        let mut stats = Stats::default();
        for &l in lens {
            stats.add_read(l, 0);
        }
        stats.n50()
    }

    #[test]
    fn n50_matches_naive() {
        let m = MAX_COUNTED_LEN;
        let cases: &[&[i64]] = &[
            &[],
            &[100],
            &[14807, 16920, 17508],          // the bam01 fixture lengths
            &[10, 10, 10, 10],               // threshold inside one bin
            &[m + 5000],                     // a single long read
            &[m + 5000, m + 5000, 10000],    // N50 among the long reads
            &[m + 5000, 50000, 50000, 50000], // N50 crosses into the short array
            &[m, m + 1],                     // boundary lengths
        ];
        for lens in cases {
            assert_eq!(stats_n50(lens), naive_n50(lens), "lens={lens:?}");
        }
    }
}

/// Whether a hit passes the alignment-length filter (`-a`). Note `-q` is NOT
/// applied here — it is a post-filter on the emitted line's col-7 mapq (see
/// `emit_trim`/`emit_split`/`emit_gap`), so it never drops an alignment upfront.
pub(crate) fn passes_filter(h: &Hit, opts: &ExtractOpts) -> bool {
    h.alen >= opts.min_aln_len
}

impl Hit {
    /// Reference offset at the read-start side of the hit.
    fn fts(&self) -> i64 {
        match self.strand {
            Strand::Fwd => self.ts,
            Strand::Rev => self.te,
        }
    }

    /// Reference offset at the read-end side of the hit.
    fn fte(&self) -> i64 {
        match self.strand {
            Strand::Fwd => self.te,
            Strand::Rev => self.ts,
        }
    }

    /// Orientation marker implied by the strand.
    fn ori(&self) -> char {
        match self.strand {
            Strand::Fwd => '>',
            Strand::Rev => '<',
        }
    }
}

fn flip(ori: char) -> char {
    match ori {
        '>' => '<',
        '<' => '>',
        c => c,
    }
}

/// A read's sequence (alignment-file input only), oriented read-forward: on the
/// original read strand, i.e. reverse-complemented (qualities reversed, not
/// complemented) when the record is reverse-mapped. It covers read-forward
/// offsets `|off, off+seq.len()|` — the whole read for a soft-clipped record, or
/// only the aligned part for a hard-clipped one (e.g. a supplementary record).
pub(crate) struct SeqView {
    pub seq: Vec<u8>,
    /// Per-base phred qualities, same orientation and length as `seq`; `None`
    /// when QUAL is missing (`*`).
    pub qual: Option<Vec<u8>>,
    pub off: i64,
}

/// Build the `;elen{sfx}=..[;equal{sfx}=..;eseq{sfx}=..]` INFO suffix for a
/// breakend window; `sfx` is `""`, or `"2"` for an insertion's second breakpoint.
///
/// The window spans read-forward query offsets `[lo, hi]` (a single point,
/// `qdist == 0`, except for a split); it is padded by `opts.flank` on each side,
/// clamped to the bases `view` covers. `elen` is emitted whenever a read sequence
/// is available (alignment input); `eseq` (on the original read strand) is
/// appended only when the window is at most `opts.max_eseq` long. `equal` (the
/// phred quality of that window) precedes `eseq` when per-base qualities are
/// available.
///
/// Returns the suffix (empty when `view` is `None`, e.g. PAF) and the window's
/// `equal`, if computed, for the caller's `-Q` check ([`below_min_equal`]).
fn eseq_info(
    view: Option<&SeqView>,
    opts: &ExtractOpts,
    lo: i64,
    hi: i64,
    qdist: i64,
    sfx: &str,
) -> (String, Option<i64>) {
    let Some(v) = view else {
        return (String::new(), None);
    };
    let start = (lo - opts.flank).max(v.off);
    let end = (hi + opts.flank).min(v.off + v.seq.len() as i64);
    let left = lo - start;
    let right = end - hi;
    let mut info = format!(";elen{sfx}={left},{qdist},{right}");
    let mut equal = None;
    if end - start <= opts.max_eseq && start <= end {
        let (s, e) = ((start - v.off) as usize, (end - v.off) as usize);
        // `equal` (phred): mean expected error rate over the window, from the
        // per-base qualities, scaled back to phred. Only when qualities exist.
        if let Some(qual) = &v.qual
            && e > s
        {
            let q = eseq_qual(&qual[s..e]);
            info.push_str(&format!(";equal{sfx}={q}"));
            equal = Some(q);
        }
        info.push_str(&format!(";eseq{sfx}="));
        info.push_str(&String::from_utf8_lossy(&v.seq[s..e]));
    }
    (info, equal)
}

/// `-Q`: whether a window's `equal` is below `opts.min_equal` (default 0 = keep
/// everything), dropping its breakend. Windows without an `equal` value pass.
fn below_min_equal(equal: Option<i64>, opts: &ExtractOpts) -> bool {
    equal.is_some_and(|q| q < opts.min_equal)
}

/// Phred quality of an eseq window: turn each base quality into an error rate,
/// average the expected errors over the window, and scale back to phred —
/// `-10 * log10( (Σ 10^(-q/10)) / L )`, rounded. For constant `q` this is `q`.
fn eseq_qual(quals: &[u8]) -> i64 {
    let expected_errors: f64 = quals
        .iter()
        .map(|&q| 10f64.powf(-(q as f64) / 10.0))
        .sum();
    let mean = expected_errors / quals.len() as f64;
    (-10.0 * mean.log10()).round() as i64
}

/// Run the extract command end to end.
pub fn run(opts: &ExtractOpts) -> io::Result<()> {
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    let mut stats = Stats::default();
    if opts.paf {
        let reader = open_reader(&opts.input)?;
        run_paf(reader, opts, &mut stats, &mut out)?;
    } else {
        crate::aln::run(&opts.input, opts, &mut stats, &mut out)?;
    }
    out.flush()?;
    stats.report();
    Ok(())
}

/// Stream a PAF `reader`, grouping hits by read name (input is assumed grouped),
/// and emit breakends for each read.
fn run_paf(
    reader: Box<dyn BufRead>,
    opts: &ExtractOpts,
    stats: &mut Stats,
    out: &mut impl Write,
) -> io::Result<()> {
    let mut group: Vec<Hit> = Vec::new();
    // Per-read stats state, tracked before the -a filter (so a fully filtered
    // read still counts): (qname, qlen, longest query span seen). PAF has no
    // primary flag, so the longest span stands in for the primary alignment.
    let mut cur: Option<(String, i64, i64)> = None;
    for line in reader.lines() {
        let line = line?;
        // An unmapped read (target `*`) has a single line and no hits; it only
        // feeds the stats.
        if let Some(qlen) = parse_paf_unmapped(&line) {
            if let Some((_, qlen, span)) = cur.take() {
                stats.add_read(qlen, span);
            }
            stats.add_unmapped(qlen);
            continue;
        }
        let Some(hit) = parse_paf_line(&line, opts.min_gap) else {
            continue;
        };
        match &mut cur {
            Some((name, _, span)) if *name == hit.qname => {
                *span = (*span).max(hit.qe - hit.qs);
            }
            _ => {
                if let Some((_, qlen, span)) = cur.take() {
                    stats.add_read(qlen, span);
                }
                cur = Some((hit.qname.clone(), hit.qlen, hit.qe - hit.qs));
            }
        }
        if !passes_filter(&hit, opts) {
            continue;
        }
        if let Some(first) = group.first()
            && first.qname != hit.qname
        {
            // PAF carries no read sequence, so eseq/elen/equal are unavailable.
            emit_read(&mut group, opts, None, out)?;
            group.clear();
        }
        group.push(hit);
    }
    emit_read(&mut group, opts, None, out)?;
    if let Some((_, qlen, span)) = cur {
        stats.add_read(qlen, span);
    }
    Ok(())
}

/// Emit the breakends of one read's hits (sorted here by `qs`): left trim (with
/// `-T`), splits (read order), right trim (with `-T`), then each hit's del/ins
/// breakends (read order; only hits whose full CIGAR is known carry gaps).
///
/// `view`, when `Some`, is the read sequence (alignment input only), enabling
/// the `elen`/`equal`/`eseq` tags. Trim/split windows need it to cover the whole
/// read; del/ins windows are clamped to what it covers.
pub(crate) fn emit_read(
    hits: &mut [Hit],
    opts: &ExtractOpts,
    view: Option<&SeqView>,
    out: &mut impl Write,
) -> io::Result<()> {
    if hits.is_empty() {
        return Ok(());
    }
    // Number of alignments for this read that passed filtering (primary +
    // supplementary). A read-level constant, stamped on every breakend it emits.
    let n_aln = hits.len() as i64;
    // Sort chimeric hits by their start position along the read.
    hits.sort_by_key(|h| h.qs);

    let first = &hits[0];
    let last = &hits[hits.len() - 1];
    // A trim/split window may reach into any part of the read.
    let full = view.filter(|v| v.off == 0 && v.seq.len() as i64 == first.qlen);

    // `idx` numbers each trim/split breakend within this read (0-based, in
    // emission order: left trim, splits, right trim); it resets for the next
    // read. A trim takes its number even when hidden (no -T), so a breakend's
    // idx does not depend on -T.
    let mut idx = 0i64;

    // Left trim: the read start is clipped and connects to nothing. Its eseq
    // window is centered on the clip point (the read start of the first hit).
    if first.qs > opts.min_clip {
        if opts.trim {
            let arrow = flip(first.ori());
            emit_trim(
                out,
                idx,
                n_aln,
                first,
                first.fts(),
                arrow,
                first.qs,
                first.qs,
                opts,
                full,
            )?;
        }
        idx += 1;
    }

    // Splits between adjacent hits along the read.
    for pair in hits.windows(2) {
        emit_split(out, idx, n_aln, &pair[0], &pair[1], opts, full)?;
        idx += 1;
    }

    // Right trim: the read end is clipped and connects to nothing. Its eseq
    // window is centered on the clip point (the read end of the last hit).
    if opts.trim && last.qlen - last.qe > opts.min_clip {
        let arrow = last.ori();
        emit_trim(
            out,
            idx,
            n_aln,
            last,
            last.fte(),
            arrow,
            last.qlen - last.qe,
            last.qe,
            opts,
            full,
        )?;
    }

    // Long indels inside each hit's CIGAR.
    for h in hits.iter() {
        emit_gaps(out, h, n_aln, view, opts)?;
    }
    Ok(())
}

/// Emit a trim breakend: a dangling end with no mate (`.` columns).
///
/// `clipped` is the unmapped read length on the dangling side; it becomes the
/// `right` field of `qlen` (the mate side), with `middle=0` and
/// `left = qlen - clipped` for the mapped hit at `ctg1:pos1`. `qpos` is the clip
/// point on the read (`qs` for a left trim, `qe` for a right trim), the center
/// of the eseq window.
#[allow(clippy::too_many_arguments)]
fn emit_trim(
    out: &mut impl Write,
    idx: i64,
    n_aln: i64,
    h: &Hit,
    pos: i64,
    arrow: char,
    clipped: i64,
    qpos: i64,
    opts: &ExtractOpts,
    view: Option<&SeqView>,
) -> io::Result<()> {
    // -q post-filter: drop the line when its col-7 mapq (the trim's mapq) is below
    // the threshold. The alignment itself is not dropped upfront.
    if h.mapq < opts.min_mapq {
        return Ok(());
    }
    let (eseq, equal) = eseq_info(view, opts, qpos, qpos, 0, "");
    if below_min_equal(equal, opts) {
        return Ok(()); // dropped by -Q (low eseq quality)
    }
    writeln!(
        out,
        "{}\t{}\t{}.\t.\t.\t{}\t{}\t{}\tsource={};ctype=trim;idx={};n_aln={};aln_len={},0;qlen={},0,{};mapq={},0{}",
        h.ctg,
        pos,
        arrow,
        h.qname,
        h.mapq,
        h.strand.as_char(),
        opts.source,
        idx,
        n_aln,
        h.alen,
        h.qlen - clipped,
        clipped,
        h.mapq,
        eseq,
    )
}

/// Emit a split breakend between upstream hit `y0` and downstream hit `y1`.
#[allow(clippy::too_many_arguments)]
fn emit_split(
    out: &mut impl Write,
    idx: i64,
    n_aln: i64,
    y0: &Hit,
    y1: &Hit,
    opts: &ExtractOpts,
    view: Option<&SeqView>,
) -> io::Result<()> {
    // Endpoints: read-end side of the upstream hit joins the read-start side
    // of the downstream hit.
    let (mut ctg0, mut pos0, mut o0, mut len0, mut mq0) =
        (y0.ctg.as_str(), y0.fte(), y0.ori(), y0.alen, y0.mapq);
    let (mut ctg1, mut pos1, mut o1, mut len1, mut mq1) =
        (y1.ctg.as_str(), y1.fts(), y1.ori(), y1.alen, y1.mapq);

    // Query lengths flanking the junction, in read order: `left` up to the end
    // of the upstream hit, `right` from the start of the downstream hit, and
    // `mid` the (possibly negative) gap between them.
    let (mut left, mut right) = (y0.qe, y0.qlen - y1.qs);
    let mid = y1.qs - y0.qe;

    // Canonicalize so the smaller (contig, pos) comes first. When we flip, both
    // orientation markers invert and the record's strand column becomes '-';
    // the alignment lengths, per-side mapqs, and the qlen left/right swap too so
    // they follow the output order (index 0 tracks ctg1:pos1, index 1 ctg2:pos2).
    let mut strand = '+';
    if !(ctg0 < ctg1 || (ctg0 == ctg1 && pos0 < pos1)) {
        let (no0, no1) = (flip(o1), flip(o0));
        std::mem::swap(&mut ctg0, &mut ctg1);
        std::mem::swap(&mut pos0, &mut pos1);
        std::mem::swap(&mut len0, &mut len1);
        std::mem::swap(&mut mq0, &mut mq1);
        std::mem::swap(&mut left, &mut right);
        o0 = no0;
        o1 = no1;
        strand = '-';
    }

    let mapq = y0.mapq.max(y1.mapq);
    // -q post-filter: drop the line when its col-7 mapq (the split's max) is below
    // the threshold. Neither alignment is dropped upfront.
    if mapq < opts.min_mapq {
        return Ok(());
    }
    // The eseq window is read-forward (not affected by the output flip): it spans
    // the junction gap/overlap between the two hits, padded by `flank`.
    let lo = y0.qe.min(y1.qs);
    let hi = y0.qe.max(y1.qs);
    let (eseq, equal) = eseq_info(view, opts, lo, hi, mid, "");
    if below_min_equal(equal, opts) {
        return Ok(()); // dropped by -Q (low eseq quality)
    }
    writeln!(
        out,
        "{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\tsource={};ctype=split;idx={};n_aln={};aln_len={},{};qlen={},{},{};mapq={},{}{}",
        ctg0,
        pos0,
        o0,
        o1,
        ctg1,
        pos1,
        y0.qname,
        mapq,
        strand,
        opts.source,
        idx,
        n_aln,
        len0,
        len1,
        left,
        mid,
        right,
        mq0,
        mq1,
        eseq,
    )
}

/// Emit the del/ins breakends of hit `h` in read-forward order; `n_aln` is the
/// read's alignment count. Called per hit by [`emit_read`], and directly for a
/// supplementary record, which contributes only its own gaps.
pub(crate) fn emit_gaps(
    out: &mut impl Write,
    h: &Hit,
    n_aln: i64,
    view: Option<&SeqView>,
    opts: &ExtractOpts,
) -> io::Result<()> {
    // The CIGAR runs along the reference, i.e. backwards along a reverse read.
    let mut gaps: Vec<&Gap> = h.gaps.iter().collect();
    if h.strand == Strand::Rev {
        gaps.reverse();
    }
    for g in gaps {
        emit_gap(out, h, g, n_aln, view, opts)?;
    }
    Ok(())
}

/// Emit one del/ins breakend from hit `h`'s CIGAR.
///
/// It is written like the split an aligner would report for the same event:
/// `ori` `>>`, the alignment strand in col 8 (what canonicalizing that split
/// gives), and `h` standing in for both hits (`aln_len=A,A`, `mapq=M,M`, col-7
/// `M`). With `q` the event's query offset in alignment orientation, `qlen` is
/// `q,0,qlen-q` (del) or `q,len,qlen-q-len` (ins) in output order.
///
/// `idx` is the event's read-forward query offset: a read's gaps come from
/// several records (primary and supplementary, in any file order), so they
/// cannot share the trim/split counter, but the offset is unique per read.
/// - del: `ctg t >> ctg t+len`; one eseq window at the read-forward junction.
/// - ins: `ctg L >> ctg R` (leftmost/rightmost equivalent placement), `inslen`,
///   and eseq windows at its two read-forward breakpoints `qa` (`elen`/`equal`/
///   `eseq`) and `qa+len` (`elen2`/`equal2`/`eseq2`), on the CIGAR placement.
fn emit_gap(
    out: &mut impl Write,
    h: &Hit,
    g: &Gap,
    n_aln: i64,
    view: Option<&SeqView>,
    opts: &ExtractOpts,
) -> io::Result<()> {
    // -q post-filter on the col-7 mapq (the hit's own).
    if h.mapq < opts.min_mapq {
        return Ok(());
    }
    let fwd = h.strand == Strand::Fwd;
    let (q, len) = (g.qoff, g.len);
    let (ctype, pos1, pos2, idx, mid, tags) = if g.ins {
        let qa = if fwd { q } else { h.qlen - q - len };
        let (e1, q1) = eseq_info(view, opts, qa, qa, 0, "");
        let (e2, q2) = eseq_info(view, opts, qa + len, qa + len, 0, "2");
        if below_min_equal(q1, opts) || below_min_equal(q2, opts) {
            return Ok(()); // dropped by -Q (low eseq quality)
        }
        ("ins", g.tpos_l, g.tpos_r, qa, len, format!(";inslen={len}{e1}{e2}"))
    } else {
        let qf = if fwd { q } else { h.qlen - q };
        let (e, q1) = eseq_info(view, opts, qf, qf, 0, "");
        if below_min_equal(q1, opts) {
            return Ok(()); // dropped by -Q (low eseq quality)
        }
        ("del", g.tpos, g.tpos + len, qf, 0, e)
    };
    writeln!(
        out,
        "{}\t{}\t>>\t{}\t{}\t{}\t{}\t{}\tsource={};ctype={};idx={};n_aln={};aln_len={},{};qlen={},{},{};mapq={},{}{}",
        h.ctg,
        pos1,
        h.ctg,
        pos2,
        h.qname,
        h.mapq,
        h.strand.as_char(),
        opts.source,
        ctype,
        idx,
        n_aln,
        h.alen,
        h.alen,
        q,
        mid,
        h.qlen - q - mid,
        h.mapq,
        h.mapq,
        tags,
    )
}
