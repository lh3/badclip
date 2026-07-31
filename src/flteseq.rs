//! The `flteseq` subcommand: filter `extract` breakends by pangenome alignment.
//!
//! Given `extract` output and the ropebwt3 `sw` PAF produced by aligning the
//! `geteseq` FASTA against a pangenome, drop each breakend whose junction is
//! already explained by the pangenome. A line is dropped if it lacks an `eseq`
//! tag, or if some ropebwt3 alignment's query interval spans (contains) the
//! junction interval on the eseq (extended by `-l` on each side) — meaning the
//! breakend is "protected" and not novel.
//!
//! The PAF's query names are `geteseq`'s names (`readName_idx_leftFlank_
//! rightFlank`) in the same order as the `eseq` lines, and a subset of them, so
//! both files are streamed with a one-line lookahead over the PAF.
//!
//! With `-s STR` (`source: Some`), the mode flips: every input line is printed,
//! and the survivors (kept/novel breakends) get their `source=` INFO tag
//! rewritten to `STR`; dropped lines are printed verbatim with their original
//! source. This tags novel calls distinctly for a downstream `merge`.

use std::io::{self, BufRead, Write};

use crate::io::open_reader;

/// Filter `extract_out` using the ropebwt3 PAF `rb3_paf`; `l` is the margin.
///
/// With `source = Some(str)`, print all input lines and rewrite the `source=`
/// tag to `str` on survivors; with `None`, print only survivors, verbatim.
pub fn run(extract_out: &str, rb3_paf: &str, l: i64, source: Option<&str>) -> io::Result<()> {
    let clip = open_reader(extract_out)?;
    let mut paf = open_reader(rb3_paf)?.lines();
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    // One-line lookahead into the PAF stream.
    let mut ahead: Option<String> = None;

    for line in clip.lines() {
        let line = line?;
        let Some((name, lo, hi)) = protected_interval(&line, l) else {
            // No eseq tag -> dropped. In print-all mode it is still emitted
            // verbatim; it has no PAF entry, so the stream stays in sync.
            if source.is_some() {
                writeln!(out, "{line}")?;
            }
            continue;
        };

        // Consume the PAF lines belonging to this eseq (same order, subset), and
        // test whether any spans [lo, hi].
        let mut protected = false;
        loop {
            if ahead.is_none() {
                ahead = match paf.next() {
                    Some(p) => Some(p?),
                    None => None,
                };
            }
            let matches = ahead
                .as_deref()
                .and_then(|p| p.split('\t').next())
                .is_some_and(|qn| qn == name);
            if !matches {
                break;
            }
            let p = ahead.take().unwrap();
            if let Some((qs, qe)) = paf_query_interval(&p)
                && qs <= lo
                && qe >= hi
            {
                protected = true;
            }
        }

        let survivor = !protected;
        match source {
            // Print-all mode: survivors get their source relabelled, dropped
            // (protected) lines are printed verbatim.
            Some(src) if survivor => writeln!(out, "{}", relabel_source(&line, src))?,
            Some(_) => writeln!(out, "{line}")?,
            // Default mode: only survivors, verbatim.
            None if survivor => writeln!(out, "{line}")?,
            None => {}
        }
    }
    Ok(())
}

/// Rewrite the `source=` value in a record's INFO (last TAB field) to `src`.
/// If no `source=` tag is present, prepend one (extract's "source first"
/// convention). Only this one tag's value changes; tag order is not significant.
fn relabel_source(line: &str, src: &str) -> String {
    let Some((rest, info)) = line.rsplit_once('\t') else {
        return line.to_string();
    };
    let mut found = false;
    let mut tags: Vec<String> = info
        .split(';')
        .map(|kv| {
            if kv.starts_with("source=") {
                found = true;
                format!("source={src}")
            } else {
                kv.to_string()
            }
        })
        .collect();
    if !found {
        tags.insert(0, format!("source={src}"));
    }
    format!("{rest}\t{}", tags.join(";"))
}

/// For an `extract` line with an `eseq` tag, return `(fastaName, lo, hi)` — the
/// geteseq name and the protected junction interval `[lo, hi]` on the eseq.
fn protected_interval(line: &str, l: i64) -> Option<(String, i64, i64)> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }
    let qname = fields[5];

    let mut idx = None;
    let mut elen = None;
    let mut has_eseq = false;
    for kv in fields[8].split(';') {
        if let Some(v) = kv.strip_prefix("idx=") {
            idx = Some(v);
        } else if let Some(v) = kv.strip_prefix("elen=") {
            elen = Some(v);
        } else if kv.starts_with("eseq=") {
            has_eseq = true;
        }
    }
    if !has_eseq {
        return None;
    }

    let (idx, elen) = (idx?, elen?);
    let e: Vec<&str> = elen.split(',').collect();
    if e.len() < 3 {
        return None;
    }
    let e0: i64 = e[0].parse().ok()?;
    let e1: i64 = e[1].parse().ok()?;
    let e2: i64 = e[2].parse().ok()?;

    let eseq_len = e0 + e1.abs() + e2;
    let lo = (e0 - l).max(0);
    let hi = (e0 + e1.abs() + l).min(eseq_len);
    // Name matches geteseq: readName_idx_leftFlank_rightFlank.
    Some((format!("{qname}_{idx}_{e0}_{e2}"), lo, hi))
}

/// Parse a PAF line's 0-based query interval `(qs, qe)` (columns 3 and 4).
fn paf_query_interval(line: &str) -> Option<(i64, i64)> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 4 {
        return None;
    }
    Some((f[2].parse().ok()?, f[3].parse().ok()?))
}
