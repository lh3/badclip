//! The `extract` subcommand: extract breakends from read alignments.
//!
//! A read may map as several chimeric hits. Sorted along the read, each read
//! end that is soft-clipped away from an alignment produces a *clip* breakend
//! (connected to nothing), and each junction between two adjacent hits
//! produces a *join* breakend. All records share one 9-column, TAB-delimited
//! layout:
//!
//! ```text
//! ctg1  pos1  ori  ctg2  pos2  qname  mapq  strand  aln_len=len1,len2
//! ```
//!
//! where `ori` is two characters, each `>`/`<`, or `.` for a missing mate. The
//! final INFO column carries `aln_len=len1,len2`, the alignment block lengths
//! (PAF col 11) of the hits at `ctg1:pos1` and `ctg2:pos2`; for a clip `len2`
//! is `0`.

use std::io::{self, BufRead, Write};

use crate::io::open_reader;
use crate::paf::{Hit, Strand, parse_paf_line};

/// Options for `badclip extract`.
pub struct ExtractOpts {
    /// Input path; `None` or `"-"` means stdin.
    pub input: Option<String>,
    /// Minimum clip length to report a clip breakend.
    pub min_clip: i64,
    /// Drop hits with mapq below this value (0 = keep everything).
    pub min_mapq: i64,
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
    let reader = open_reader(opts.input.as_deref())?;
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    extract(reader, opts, &mut out)
}

/// Stream `reader`, grouping hits by read name (input is assumed grouped),
/// and emit breakends for each read.
fn extract(reader: Box<dyn BufRead>, opts: &ExtractOpts, out: &mut impl Write) -> io::Result<()> {
    let mut group: Vec<Hit> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let Some(hit) = parse_paf_line(&line) else {
            continue;
        };
        if hit.mapq < opts.min_mapq {
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

/// Emit clip and join breakends for one read's hits.
fn emit_read(hits: &mut [Hit], opts: &ExtractOpts, out: &mut impl Write) -> io::Result<()> {
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
        emit_clip(out, first, first.fts(), arrow)?;
    }

    // Joins between adjacent hits along the read.
    for pair in hits.windows(2) {
        emit_join(out, &pair[0], &pair[1])?;
    }

    // Right clip: the read end is clipped and connects to nothing.
    if last.qlen - last.qe > opts.min_clip {
        let arrow = last.ori();
        emit_clip(out, last, last.fte(), arrow)?;
    }

    Ok(())
}

/// Emit a clip breakend: a dangling end with no mate (`.` columns).
fn emit_clip(out: &mut impl Write, h: &Hit, pos: i64, arrow: char) -> io::Result<()> {
    writeln!(
        out,
        "{}\t{}\t{}.\t.\t.\t{}\t{}\t{}\taln_len={},0",
        h.ctg,
        pos,
        arrow,
        h.qname,
        h.mapq,
        h.strand.as_char(),
        h.alen,
    )
}

/// Emit a join breakend between upstream hit `y0` and downstream hit `y1`.
fn emit_join(out: &mut impl Write, y0: &Hit, y1: &Hit) -> io::Result<()> {
    // Endpoints: read-end side of the upstream hit joins the read-start side
    // of the downstream hit.
    let (mut ctg0, mut pos0, mut o0, mut len0) = (y0.ctg.as_str(), y0.fte(), y0.ori(), y0.alen);
    let (mut ctg1, mut pos1, mut o1, mut len1) = (y1.ctg.as_str(), y1.fts(), y1.ori(), y1.alen);

    // Canonicalize so the smaller (contig, pos) comes first. When we flip, both
    // orientation markers invert and the record's strand column becomes '-';
    // the alignment lengths swap too so `aln_len` follows the output order.
    let mut strand = '+';
    if !(ctg0 < ctg1 || (ctg0 == ctg1 && pos0 < pos1)) {
        let (no0, no1) = (flip(o1), flip(o0));
        std::mem::swap(&mut ctg0, &mut ctg1);
        std::mem::swap(&mut pos0, &mut pos1);
        std::mem::swap(&mut len0, &mut len1);
        o0 = no0;
        o1 = no1;
        strand = '-';
    }

    let mapq = y0.mapq.min(y1.mapq);
    writeln!(
        out,
        "{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\taln_len={},{}",
        ctg0, pos0, o0, o1, ctg1, pos1, y0.qname, mapq, strand, len0, len1,
    )
}
