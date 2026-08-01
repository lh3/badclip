//! Alignment-file input for `extract` (BAM, CRAM, and SAM).
//!
//! Format is autodetected by noodles' `noodles_util` alignment reader
//! (SAM/BAM/CRAM by magic bytes), so one path serves all three. A read's
//! chimeric hits are recovered from the primary alignment alone: its `SA:Z:` tag
//! lists the supplementary alignments. So we iterate primary records only —
//! skipping secondary, supplementary, and unmapped records — and reconstruct
//! every hit from the primary's CIGAR plus the `SA` entries. This needs no record
//! grouping, so it works on both name-grouped and coordinate-sorted files.
//!
//! CRAM stores sequences as diffs against a reference, so it requires the
//! reference genome (`-r`, a faidx-indexed FASTA); the reader decodes each
//! record's sequence against it. BAM/SAM carry the sequence directly and ignore
//! `-r`.
//!
//! Hit coordinates are computed to match the PAF path exactly (see
//! [`span_from_ops`]). The only value that can differ from PAF is `aln_len` for
//! supplementary hits, because minimap2's `SA` CIGAR is collapsed (it merges
//! small indels); the primary hit uses its full CIGAR and matches PAF.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};

use noodles_fasta as fasta;
use noodles_sam::Header;
use noodles_sam::alignment::Record;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_util::alignment;

use crate::extract::{ExtractOpts, emit_read, passes_filter};
use crate::paf::{Hit, Strand};

/// Read an alignment file from `input` (`"-"` = stdin), autodetecting
/// SAM/BAM/CRAM, and emit breakends. CRAM requires `opts.reference` (`-r`).
pub fn run(input: &str, opts: &ExtractOpts, out: &mut impl Write) -> io::Result<()> {
    let inner: Box<dyn Read> = if input == "-" {
        Box::new(io::stdin().lock())
    } else {
        Box::new(File::open(input)?)
    };
    // Peek the magic bytes (non-consuming) only to gate the reference
    // requirement; noodles-util still performs the real format detection.
    let mut reader = BufReader::new(inner);
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

    for result in reader.records(&header) {
        emit_record(&*result?, &header, &names, opts, out)?;
    }
    Ok(())
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

/// Turn one alignment record into breakends via the shared [`emit_read`].
fn emit_record(
    record: &dyn Record,
    header: &Header,
    names: &[String],
    opts: &ExtractOpts,
    out: &mut impl Write,
) -> io::Result<()> {
    let flags = record.flags()?;
    if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped() {
        return Ok(());
    }

    let Some(ref_id) = record.reference_sequence_id(header).transpose()? else {
        return Ok(());
    };
    let Some(ctg) = names.get(ref_id) else {
        return Ok(());
    };
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

    // Primary hit, from the record's own (detailed) CIGAR.
    let mut ops: Vec<(char, i64)> = Vec::new();
    for op in record.cigar().iter() {
        let op = op?;
        ops.push((kind_char(op.kind()), op.len() as i64));
    }
    let primary = make_hit(qname.clone(), ctg.clone(), pos0, strand, mapq, &ops);
    let qlen = primary.qlen;
    let mut hits = vec![primary];

    // Supplementary hits, from the SA:Z: tag.
    if let Some(sa) = sa_tag(record)? {
        for entry in sa.split(';') {
            if let Some(hit) = parse_sa_entry(&qname, entry) {
                hits.push(hit);
            }
        }
    }

    // Read sequence on the original read strand, for eseq extraction. The record
    // SEQ is in alignment orientation (decoded against the reference for CRAM);
    // reverse-complement it when the primary is reverse-mapped. Only usable when
    // it spans the whole read (i.e. the primary is soft-, not hard-, clipped).
    let seq: Vec<u8> = record.sequence().iter().collect();
    let read_fwd: Option<Vec<u8>> = if seq.len() as i64 == qlen {
        Some(match strand {
            Strand::Rev => revcomp(&seq),
            Strand::Fwd => seq,
        })
    } else {
        None
    };

    // Per-base qualities on the original read strand, for the `equal` tag. Raw
    // phred values; on the reverse strand they are reversed (NOT complemented) to
    // stay aligned with `read_fwd`. Absent when QUAL is missing (`*`).
    let qual: Vec<u8> = record
        .quality_scores()
        .iter()
        .collect::<io::Result<Vec<u8>>>()?;
    let read_fwd_qual: Option<Vec<u8>> = if qual.len() as i64 == qlen {
        Some(match strand {
            Strand::Rev => qual.into_iter().rev().collect(),
            Strand::Fwd => qual,
        })
    } else {
        None
    };

    hits.retain(|h| passes_filter(h, opts));
    emit_read(
        &mut hits,
        opts,
        read_fwd.as_deref(),
        read_fwd_qual.as_deref(),
        out,
    )
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
    Some(make_hit(qname.to_string(), ctg, pos0, strand, mapq, &ops))
}

/// Build a [`Hit`] from a hit's reference start, strand, mapq and CIGAR ops.
fn make_hit(
    qname: String,
    ctg: String,
    pos0: i64,
    strand: Strand,
    mapq: i64,
    ops: &[(char, i64)],
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

/// Parse a CIGAR string into `(op_char, len)` pairs.
fn parse_cigar(s: &str) -> Vec<(char, i64)> {
    let mut ops = Vec::new();
    let mut num: i64 = 0;
    for b in s.bytes() {
        if b.is_ascii_digit() {
            num = num * 10 + (b - b'0') as i64;
        } else {
            ops.push((b as char, num));
            num = 0;
        }
    }
    ops
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
