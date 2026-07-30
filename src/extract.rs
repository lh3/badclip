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
            emit_read(&mut group, opts, out)?;
            group.clear();
        }
        group.push(hit);
    }
    emit_read(&mut group, opts, out)?;
    Ok(())
}

/// Emit clip and join breakends for one read's hits (sorted here by `qs`).
pub(crate) fn emit_read(
    hits: &mut [Hit],
    opts: &ExtractOpts,
    out: &mut impl Write,
) -> io::Result<()> {
    if hits.is_empty() {
        return Ok(());
    }
    // Sort chimeric hits by their start position along the read.
    hits.sort_by_key(|h| h.qs);

    let first = &hits[0];
    let last = &hits[hits.len() - 1];

    // Left clip: the read start is clipped and connects to nothing.
    if first.qs > opts.min_clip {
        let arrow = flip(first.ori());
        emit_clip(out, first, first.fts(), arrow, first.qs)?;
    }

    // Joins between adjacent hits along the read.
    for pair in hits.windows(2) {
        emit_join(out, &pair[0], &pair[1])?;
    }

    // Right clip: the read end is clipped and connects to nothing.
    if last.qlen - last.qe > opts.min_clip {
        let arrow = last.ori();
        emit_clip(out, last, last.fte(), arrow, last.qlen - last.qe)?;
    }

    Ok(())
}

/// Emit a clip breakend: a dangling end with no mate (`.` columns).
///
/// `clipped` is the unmapped read length on the dangling side; it becomes the
/// `right` field of `qlen` (the mate side), with `middle=0` and
/// `left = qlen - clipped` for the mapped hit at `ctg1:pos1`.
fn emit_clip(out: &mut impl Write, h: &Hit, pos: i64, arrow: char, clipped: i64) -> io::Result<()> {
    writeln!(
        out,
        "{}\t{}\t{}.\t.\t.\t{}\t{}\t{}\taln_len={},0;qlen={},0,{};mapq={},0",
        h.ctg,
        pos,
        arrow,
        h.qname,
        h.mapq,
        h.strand.as_char(),
        h.alen,
        h.qlen - clipped,
        clipped,
        h.mapq,
    )
}

/// Emit a join breakend between upstream hit `y0` and downstream hit `y1`.
fn emit_join(out: &mut impl Write, y0: &Hit, y1: &Hit) -> io::Result<()> {
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
    writeln!(
        out,
        "{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\taln_len={},{};qlen={},{},{};mapq={},{}",
        ctg0,
        pos0,
        o0,
        o1,
        ctg1,
        pos1,
        y0.qname,
        mapq,
        strand,
        len0,
        len1,
        left,
        mid,
        right,
        mq0,
        mq1,
    )
}
