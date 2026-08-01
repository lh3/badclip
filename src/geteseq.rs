//! The `geteseq` subcommand: turn `extract` output into a FASTA of the
//! extracted breakend sequences.
//!
//! Reads `extract`'s TAB-delimited records and, for each one that carries an
//! `eseq` INFO tag, writes a FASTA record whose name is
//! `readName_idx_leftFlank_rightFlank` (from the `idx` and `elen` tags) and
//! whose sequence is the `eseq` value. Records without `eseq` are skipped, as are
//! records whose `equal` quality is below `-Q` (records lacking `equal` are kept).

use std::io::{self, BufRead, Write};

use crate::io::open_reader;

/// Read `extract` output from `input` (`"-"` = stdin) and emit FASTA. Records
/// whose `equal` tag is below `min_equal` are dropped.
pub fn run(input: &str, min_equal: i64) -> io::Result<()> {
    let reader = open_reader(input)?;
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    for line in reader.lines() {
        let line = line?;
        if let Some((name, seq)) = fasta_record(&line, min_equal) {
            writeln!(out, ">{name}")?;
            writeln!(out, "{seq}")?;
        }
    }
    Ok(())
}

/// Parse one `extract` record; return `(name, seq)` when it has an `eseq` tag and
/// its `equal` quality is `>= min_equal` (records without `equal` pass).
fn fasta_record(line: &str, min_equal: i64) -> Option<(String, &str)> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }
    let qname = fields[5];

    let mut idx = None;
    let mut elen = None;
    let mut eseq = None;
    let mut equal = None;
    for kv in fields[8].split(';') {
        if let Some(v) = kv.strip_prefix("idx=") {
            idx = Some(v);
        } else if let Some(v) = kv.strip_prefix("elen=") {
            elen = Some(v);
        } else if let Some(v) = kv.strip_prefix("eseq=") {
            eseq = Some(v);
        } else if let Some(v) = kv.strip_prefix("equal=") {
            equal = v.parse::<i64>().ok();
        }
    }

    // Drop records with a low eseq quality; keep those lacking an `equal` tag.
    if equal.is_some_and(|q| q < min_equal) {
        return None;
    }

    let (eseq, idx, elen) = (eseq?, idx?, elen?);
    let e: Vec<&str> = elen.split(',').collect();
    if e.len() < 3 {
        return None;
    }
    // Name: readName_idx_leftFlank_rightFlank (elen[0], elen[2]).
    Some((format!("{qname}_{idx}_{}_{}", e[0], e[2]), eseq))
}
