//! The `extract` subcommand: extract breakends from read alignments.
//!
//! A read may map as several chimeric hits. Sorted along the read, each read
//! end that is soft-clipped away from an alignment produces a *clip* breakend
//! (connected to nothing), and each junction between two adjacent hits
//! produces a *join* breakend. All records share one 9-column, TAB-delimited
//! layout:
//!
//! ```text
//! ctg1  pos1  ori  ctg2  pos2  qname  mapq  strand  aln_len=..;qlen=..;mapq=..
//! ```
//!
//! where `ori` is two characters, each `>`/`<`, or `.` for a missing mate. The
//! final INFO column carries three tags:
//! - `aln_len=len1,len2` — alignment block lengths (PAF col 11) of the hits at
//!   `ctg1:pos1` and `ctg2:pos2`; for a clip `len2` is `0`.
//! - `qlen=left,middle,right` — query lengths whose sum is the read length.
//!   `left` is on the `ctg1:pos1` side, `right` on the `ctg2:pos2` side. For a
//!   clip, `right` is the clipped length and `middle` is `0`; for a join,
//!   `middle` is the (possibly negative) query gap between the two hits.
//! - `mapq=mapq1,mapq2` — the mapq of each hit, `mapq1` for `ctg1:pos1` and
//!   `mapq2` for `ctg2:pos2`; for a clip `mapq2` is `0`. (Distinct from the
//!   field-7 mapq column, which is the `min` of a join's two mapqs.)

use std::io::{self, BufRead, Write};

use crate::io::open_reader;
use crate::paf::{Hit, Strand, parse_paf_line};

/// Options for `badclip extract`.
pub struct ExtractOpts {
    /// Input path; `"-"` means stdin.
    pub input: String,
    /// Read PAF instead of BAM.
    pub paf: bool,
    /// Minimum clip length to report a clip breakend.
    pub min_clip: i64,
    /// Drop hits with mapq below this value (0 = keep everything).
    pub min_mapq: i64,
    /// Drop hits whose alignment block length is below this (0 = keep all).
    pub min_aln_len: i64,
    /// Flanking read sequence extracted on each side of a breakend (BAM only).
    pub flank: i64,
    /// Maximum extracted sequence length; longer windows omit `eseq` (BAM only).
    pub max_eseq: i64,
}

/// Whether a hit passes the mapq / alignment-length filters.
pub(crate) fn passes_filter(h: &Hit, opts: &ExtractOpts) -> bool {
    h.mapq >= opts.min_mapq && h.alen >= opts.min_aln_len
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

/// Build the `;elen=..[;eseq=..]` INFO suffix for a breakend window.
///
/// The window spans read-forward query offsets `[lo, hi]` (a single point for a
/// clip, `qdist == 0`); it is padded by `opts.flank` on each side, clamped to
/// the read. `elen` is emitted whenever a read sequence is available (BAM);
/// `eseq` (the substring `read_seq[start..end]`, on the original read strand) is
/// appended only when the window is at most `opts.max_eseq` long. Returns `""`
/// when `read_seq` is `None` (PAF).
fn eseq_info(
    read_seq: Option<&[u8]>,
    opts: &ExtractOpts,
    qlen: i64,
    lo: i64,
    hi: i64,
    qdist: i64,
) -> String {
    let Some(seq) = read_seq else {
        return String::new();
    };
    let start = (lo - opts.flank).max(0);
    let end = (hi + opts.flank).min(qlen);
    let left = lo - start;
    let right = end - hi;
    let mut info = format!(";elen={left},{qdist},{right}");
    if end - start <= opts.max_eseq && start >= 0 && start <= end && (end as usize) <= seq.len() {
        info.push_str(";eseq=");
        info.push_str(&String::from_utf8_lossy(&seq[start as usize..end as usize]));
    }
    info
}

/// Run the extract command end to end.
pub fn run(opts: &ExtractOpts) -> io::Result<()> {
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    if opts.paf {
        let reader = open_reader(&opts.input)?;
        run_paf(reader, opts, &mut out)
    } else {
        crate::bam::run_bam(&opts.input, opts, &mut out)
    }
}

/// Stream a PAF `reader`, grouping hits by read name (input is assumed grouped),
/// and emit breakends for each read.
fn run_paf(reader: Box<dyn BufRead>, opts: &ExtractOpts, out: &mut impl Write) -> io::Result<()> {
    let mut group: Vec<Hit> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let Some(hit) = parse_paf_line(&line) else {
            continue;
        };
        if !passes_filter(&hit, opts) {
            continue;
        }
        if let Some(first) = group.first()
            && first.qname != hit.qname
        {
            // PAF carries no read sequence, so eseq/elen are unavailable.
            emit_read(&mut group, opts, None, out)?;
            group.clear();
        }
        group.push(hit);
    }
    emit_read(&mut group, opts, None, out)?;
    Ok(())
}

/// Emit clip and join breakends for one read's hits (sorted here by `qs`).
///
/// `read_seq`, when `Some`, is the read sequence on the original read strand
/// (BAM only), enabling the `elen`/`eseq` tags.
pub(crate) fn emit_read(
    hits: &mut [Hit],
    opts: &ExtractOpts,
    read_seq: Option<&[u8]>,
    out: &mut impl Write,
) -> io::Result<()> {
    if hits.is_empty() {
        return Ok(());
    }
    // Sort chimeric hits by their start position along the read.
    hits.sort_by_key(|h| h.qs);

    let first = &hits[0];
    let last = &hits[hits.len() - 1];

    // `idx` numbers each emitted breakend within this read (0-based, in emission
    // order: left clip, joins, right clip); it resets for the next read.
    let mut idx = 0i64;

    // Left clip: the read start is clipped and connects to nothing. Its eseq
    // window is centered on the clip point (the read start of the first hit).
    if first.qs > opts.min_clip {
        let arrow = flip(first.ori());
        emit_clip(
            out,
            idx,
            first,
            first.fts(),
            arrow,
            first.qs,
            first.qs,
            opts,
            read_seq,
        )?;
        idx += 1;
    }

    // Joins between adjacent hits along the read.
    for pair in hits.windows(2) {
        emit_join(out, idx, &pair[0], &pair[1], opts, read_seq)?;
        idx += 1;
    }

    // Right clip: the read end is clipped and connects to nothing. Its eseq
    // window is centered on the clip point (the read end of the last hit).
    if last.qlen - last.qe > opts.min_clip {
        let arrow = last.ori();
        emit_clip(
            out,
            idx,
            last,
            last.fte(),
            arrow,
            last.qlen - last.qe,
            last.qe,
            opts,
            read_seq,
        )?;
    }

    Ok(())
}

/// Emit a clip breakend: a dangling end with no mate (`.` columns).
///
/// `clipped` is the unmapped read length on the dangling side; it becomes the
/// `right` field of `qlen` (the mate side), with `middle=0` and
/// `left = qlen - clipped` for the mapped hit at `ctg1:pos1`. `qpos` is the clip
/// point on the read (`qs` for a left clip, `qe` for a right clip), the center
/// of the eseq window.
#[allow(clippy::too_many_arguments)]
fn emit_clip(
    out: &mut impl Write,
    idx: i64,
    h: &Hit,
    pos: i64,
    arrow: char,
    clipped: i64,
    qpos: i64,
    opts: &ExtractOpts,
    read_seq: Option<&[u8]>,
) -> io::Result<()> {
    let eseq = eseq_info(read_seq, opts, h.qlen, qpos, qpos, 0);
    writeln!(
        out,
        "{}\t{}\t{}.\t.\t.\t{}\t{}\t{}\tidx={};aln_len={},0;qlen={},0,{};mapq={},0{}",
        h.ctg,
        pos,
        arrow,
        h.qname,
        h.mapq,
        h.strand.as_char(),
        idx,
        h.alen,
        h.qlen - clipped,
        clipped,
        h.mapq,
        eseq,
    )
}

/// Emit a join breakend between upstream hit `y0` and downstream hit `y1`.
#[allow(clippy::too_many_arguments)]
fn emit_join(
    out: &mut impl Write,
    idx: i64,
    y0: &Hit,
    y1: &Hit,
    opts: &ExtractOpts,
    read_seq: Option<&[u8]>,
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

    let mapq = y0.mapq.min(y1.mapq);
    // The eseq window is read-forward (not affected by the output flip): it spans
    // the junction gap/overlap between the two hits, padded by `flank`.
    let lo = y0.qe.min(y1.qs);
    let hi = y0.qe.max(y1.qs);
    let eseq = eseq_info(read_seq, opts, y0.qlen, lo, hi, mid);
    writeln!(
        out,
        "{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\tidx={};aln_len={},{};qlen={},{},{};mapq={},{}{}",
        ctg0,
        pos0,
        o0,
        o1,
        ctg1,
        pos1,
        y0.qname,
        mapq,
        strand,
        idx,
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
