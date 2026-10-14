//! Alignment-file input for `extract` (BAM, CRAM, and SAM).
//!
//! Format is autodetected by noodles' `noodles_util` alignment reader
//! (SAM/BAM/CRAM by magic bytes), so one path serves all three. A read's
//! chimeric hits are recovered from the primary alignment alone: its `SA:Z:` tag
//! lists the supplementary alignments. So the trims/splits come from primary
//! records only — every hit reconstructed from the primary's CIGAR plus the `SA`
//! entries — while secondary records are skipped and unmapped ones only feed the
//! end-of-run stats. The `SA` CIGARs are collapsed (indels merged), so long
//! CIGAR indels (del/ins) come from each record's own CIGAR: a supplementary
//! record contributes just its own del/ins breakends. Each record is handled on
//! its own, so this needs no record grouping and works on both name-grouped and
//! coordinate-sorted files.
//!
//! CRAM stores sequences as diffs against a reference, so it requires the
//! reference genome (`-r`, a faidx-indexed FASTA); the reader decodes each
//! record's sequence against it. BAM/SAM carry the sequence directly and ignore
//! `-r`.
//!
//! Hit coordinates are computed to match the PAF path exactly (see
//! [`span_from_ops`]). Two things can differ from PAF: `aln_len` for
//! supplementary hits, because minimap2's `SA` CIGAR is collapsed (it merges
//! small indels) — and minimap2's PAF col 11 can differ from a full CIGAR's
//! M+I+D on alignments with long indels — and insertion positions, which are
//! normalized with the read sequence (absent from PAF).

use std::collections::HashSet;
use std::io::{self, BufRead, BufReader, Write};

use noodles_fasta as fasta;
use noodles_sam::Header;
use noodles_sam::alignment::Record;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_util::alignment;

use crate::extract::{ExtractOpts, SeqView, Stats, emit_gaps, emit_read, passes_filter};
use crate::io::open_reader;
use crate::paf::{Gap, Hit, Strand, find_gaps, parse_cigar};

/// Read an alignment file from `input` (`"-"` = stdin), autodetecting
/// SAM/BAM/CRAM, and emit breakends. CRAM requires `opts.reference` (`-r`).
pub fn run(
    input: &str,
    opts: &ExtractOpts,
    stats: &mut Stats,
    out: &mut impl Write,
) -> io::Result<()> {
    // Peek the magic bytes (non-consuming) only to gate the reference
    // requirement; noodles-util still performs the real format detection.
    let mut reader = BufReader::new(crate::io::open_raw(input)?);
    let is_cram = reader.fill_buf()?.starts_with(b"CRAM");

    let mut builder = alignment::io::reader::Builder::default();
    if is_cram {
        let Some(reference) = opts.reference.as_deref() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "CRAM input requires a reference genome (-r)",
            ));
        };
        builder = builder.set_reference_sequence_repository(reference_repository(reference)?);
    }

    let mut reader = builder.build_from_reader(reader)?;
    let header = reader.read_header()?;
    // Reference names, indexed by reference_sequence_id.
    let names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| k.to_string())
        .collect();

    let alt = load_alt(opts.alt.as_deref())?;

    for result in reader.records(&header) {
        emit_record(&*result?, &header, &names, &alt, opts, stats, out)?;
    }
    Ok(())
}

/// Load the set of ALT contig names from `--alt` (gzip ok). Lines starting with
/// `@` (SAM headers) are ignored; for the rest the first column (tab-delimited)
/// is taken as the contig name. So both a plain one-name-per-line list and a
/// bwa-kit `.alt` file (SAM records whose QNAME is the ALT contig, possibly with
/// `@` headers) work. Empty when no file is given.
fn load_alt(path: Option<&str>) -> io::Result<HashSet<String>> {
    let mut set = HashSet::new();
    let Some(path) = path else {
        return Ok(set);
    };
    for line in open_reader(path)?.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }
        if let Some(name) = line.split('\t').next()
            && !name.is_empty()
        {
            set.insert(name.to_string());
        }
    }
    Ok(set)
}

/// Build a CRAM reference repository from a faidx-indexed FASTA (`ref` + `ref.fai`).
fn reference_repository(reference: &str) -> io::Result<fasta::Repository> {
    let indexed = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("reading reference `{reference}` (must be faidx-indexed, i.e. `{reference}.fai` must exist): {e}"),
            )
        })?;
    Ok(fasta::Repository::new(
        fasta::repository::adapters::IndexedReader::new(indexed),
    ))
}

/// Turn one alignment record into breakends.
///
/// A **primary** record yields the read's trims/splits (from its own hit plus
/// the `SA` hits) and its own del/ins breakends, via the shared [`emit_read`].
/// A **supplementary** record yields only its own del/ins breakends: the `SA`
/// CIGARs are collapsed, so a supplementary indel is visible only here.
/// Secondary records are ignored.
///
/// `alt` is the set of ALT contig names (`--alt`): a read whose primary hit lands
/// on an ALT contig is skipped entirely (for a supplementary record, the primary
/// is the first `SA` entry, by SAM convention), and ALT hits are dropped before
/// breakends are formed.
#[allow(clippy::too_many_arguments)]
fn emit_record(
    record: &dyn Record,
    header: &Header,
    names: &[String],
    alt: &HashSet<String>,
    opts: &ExtractOpts,
    stats: &mut Stats,
    out: &mut impl Write,
) -> io::Result<()> {
    let flags = record.flags()?;
    if flags.is_secondary() {
        return Ok(());
    }
    let supp = flags.is_supplementary();
    // Unmapped primary records only feed the end-of-run stats.
    if flags.is_unmapped() {
        if !supp {
            stats.add_unmapped(record.sequence().len() as i64);
        }
        return Ok(());
    }

    let Some(ref_id) = record.reference_sequence_id(header).transpose()? else {
        return Ok(());
    };
    let Some(ctg) = names.get(ref_id) else {
        return Ok(());
    };
    // Ignore reads whose primary hit lands on an ALT contig (and supplementary
    // hits on ALT contigs).
    if alt.contains(ctg.as_str()) {
        return Ok(());
    }
    let Some(pos) = record.alignment_start().transpose()? else {
        return Ok(());
    };
    let Some(name) = record.name() else {
        return Ok(());
    };
    // For paired-end reads, disambiguate the two mates by appending `/1` (first
    // segment) or `/2` (last segment). Unpaired reads (e.g. long reads) keep the
    // bare name. The suffix flows into the SA-derived hits via `qname` too.
    let mut qname = name.to_string();
    if flags.is_segmented() {
        if flags.is_first_segment() {
            qname.push_str("/1");
        } else if flags.is_last_segment() {
            qname.push_str("/2");
        }
    }
    let pos0 = pos.get() as i64 - 1;
    let strand = strand_of(flags.is_reverse_complemented());
    let mapq = record
        .mapping_quality()
        .transpose()?
        .map_or(255, |q| q.get() as i64);

    // This record's own hit, from its (detailed) CIGAR.
    let mut ops: Vec<(char, i64)> = Vec::new();
    for op in record.cigar().iter() {
        let op = op?;
        ops.push((kind_char(op.kind()), op.len() as i64));
    }
    // A supplementary record matters only for its long indels; skip the rest
    // (incl. decoding its sequence) when it has none.
    if supp && find_gaps(&ops, 0, pos0, opts.min_gap, None).is_empty() {
        return Ok(());
    }

    // Supplementary hits, from the SA:Z: tag. Skip supplementary hits on ALT
    // contigs so they never form breakends with the (non-ALT) primary.
    let mut sa_hits = Vec::new();
    if let Some(sa) = sa_tag(record)? {
        for (i, entry) in sa.split(';').enumerate() {
            let Some(hit) = parse_sa_entry(&qname, entry) else {
                continue;
            };
            if alt.contains(hit.ctg.as_str()) {
                // A supplementary record's first SA entry is its primary.
                if supp && i == 0 {
                    return Ok(());
                }
                continue;
            }
            sa_hits.push(hit);
        }
    }

    // The stored SEQ (alignment orientation, decoded against the reference for
    // CRAM) normalizes insertions; its first base is at the leading hard clip.
    let seq: Vec<u8> = record.sequence().iter().collect();
    let lead_hard = match ops.first() {
        Some(&('H', len)) => len,
        _ => 0,
    };
    let aligned_seq = (!seq.is_empty()).then_some((seq.as_slice(), lead_hard));
    let gaps = find_gaps(&ops, 0, pos0, opts.min_gap, aligned_seq);
    let hit = make_hit(qname, ctg.clone(), pos0, strand, mapq, &ops, gaps);
    let qlen = hit.qlen;
    if !supp {
        // Count the read before any filtering; the primary's query span
        // (qe - qs) is its contribution to the primary-aligned base total.
        stats.add_read(qlen, hit.qe - hit.qs);
    }

    // Per-base qualities, raw phred values in the same orientation as SEQ.
    // Absent when QUAL is missing (`*`).
    let qual: Vec<u8> = record
        .quality_scores()
        .iter()
        .collect::<io::Result<Vec<u8>>>()?;
    let view = seq_view(seq, qual, strand, qlen, hit.qs, hit.qe);

    sa_hits.retain(|h| passes_filter(h, opts));
    if supp {
        if !passes_filter(&hit, opts) {
            return Ok(());
        }
        let n_aln = 1 + sa_hits.len() as i64;
        return emit_gaps(out, &hit, n_aln, view.as_ref(), opts);
    }
    let mut hits = vec![hit];
    hits.retain(|h| passes_filter(h, opts));
    hits.extend(sa_hits);
    emit_read(&mut hits, opts, view.as_ref(), out)
}

/// Orient a record's SEQ/QUAL read-forward (the original read strand): the
/// record SEQ is in alignment orientation, so reverse-complement it — and
/// reverse (NOT complement) the qualities — when the record is reverse-mapped.
/// SEQ must span either the whole read (soft-clipped record; offset 0) or just
/// the aligned part `|qs,qe|` (hard-clipped; offset `qs`), else `None`. The
/// qualities are kept only when they match SEQ in length.
fn seq_view(
    seq: Vec<u8>,
    qual: Vec<u8>,
    strand: Strand,
    qlen: i64,
    qs: i64,
    qe: i64,
) -> Option<SeqView> {
    let n = seq.len() as i64;
    let off = if n == qlen {
        0
    } else if n > 0 && n == qe - qs {
        qs
    } else {
        return None;
    };
    let qual = (qual.len() == seq.len()).then(|| match strand {
        Strand::Rev => qual.into_iter().rev().collect(),
        Strand::Fwd => qual,
    });
    let seq = match strand {
        Strand::Rev => revcomp(&seq),
        Strand::Fwd => seq,
    };
    Some(SeqView { seq, qual, off })
}

/// Reverse-complement a read sequence (ASCII bases).
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Complement a single ASCII base (ACGT + N + IUPAC codes; unknown pass through).
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    }
}

fn strand_of(reverse: bool) -> Strand {
    if reverse { Strand::Rev } else { Strand::Fwd }
}

/// Read the `SA:Z:` tag value as an owned string, if present.
fn sa_tag(record: &dyn Record) -> io::Result<Option<String>> {
    match record.data().get(&Tag::new(b'S', b'A')) {
        Some(value) => match value? {
            Value::String(s) => Ok(Some(s.to_string())),
            _ => Ok(None),
        },
        None => Ok(None),
    }
}

/// Parse one `SA` entry `rname,pos,strand,CIGAR,mapQ,NM` into a [`Hit`].
fn parse_sa_entry(qname: &str, entry: &str) -> Option<Hit> {
    if entry.is_empty() {
        return None;
    }
    let f: Vec<&str> = entry.split(',').collect();
    if f.len() < 6 {
        return None;
    }
    let ctg = f[0].to_string();
    let pos0 = f[1].parse::<i64>().ok()? - 1;
    let strand = strand_of(f[2] == "-");
    let ops = parse_cigar(f[3]);
    let mapq = f[4].parse::<i64>().ok()?;
    // The SA CIGAR is collapsed (indels merged), so its gaps are unknown.
    Some(make_hit(
        qname.to_string(),
        ctg,
        pos0,
        strand,
        mapq,
        &ops,
        Vec::new(),
    ))
}

/// Build a [`Hit`] from a hit's reference start, strand, mapq, CIGAR ops and
/// long indels.
fn make_hit(
    qname: String,
    ctg: String,
    pos0: i64,
    strand: Strand,
    mapq: i64,
    ops: &[(char, i64)],
    gaps: Vec<Gap>,
) -> Hit {
    let (qs, qe, ts, te, alen, qlen) = span_from_ops(pos0, strand, ops);
    Hit {
        qname,
        qlen,
        qs,
        qe,
        strand,
        ctg,
        ts,
        te,
        alen,
        mapq,
        gaps,
    }
}

/// Derive `(qs, qe, ts, te, alen, qlen)` from a CIGAR, in read-forward query
/// coordinates (identical to the PAF fields). `ops` are `(op_char, len)` pairs.
fn span_from_ops(pos0: i64, strand: Strand, ops: &[(char, i64)]) -> (i64, i64, i64, i64, i64, i64) {
    let mut qspan = 0; // query-consumed by the alignment (M/=/X/I)
    let mut refspan = 0; // reference-consumed (M/=/X/D/N)
    let mut alen = 0; // alignment block length (M/=/X/I/D)
    for &(op, len) in ops {
        match op {
            'M' | '=' | 'X' => {
                qspan += len;
                refspan += len;
                alen += len;
            }
            'I' => {
                qspan += len;
                alen += len;
            }
            'D' => {
                refspan += len;
                alen += len;
            }
            'N' => refspan += len,
            _ => {} // S/H (clips) and P (pad): not aligned
        }
    }
    let clip = |&(op, len): &(char, i64)| if op == 'S' || op == 'H' { len } else { 0 };
    let lead = ops.first().map_or(0, clip);
    let tail = ops.last().map_or(0, clip);
    let qlen = lead + qspan + tail;
    // For '-' the read is reverse-complemented, so the leading/trailing clips
    // swap when mapped back to read-forward coordinates.
    let qs = if strand == Strand::Fwd { lead } else { tail };
    (qs, qs + qspan, pos0, pos0 + refspan, alen, qlen)
}

/// Map a noodles CIGAR op kind to its SAM character.
fn kind_char(kind: Kind) -> char {
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}
