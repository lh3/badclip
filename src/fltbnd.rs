//! The `fltbnd` subcommand (alias `flteseq`): filter `extract` breakends by
//! pangenome alignment.
//!
//! Given `extract` output and the ropebwt3 `sw` PAF produced by aligning the
//! `getfa` FASTA against a pangenome, drop each breakend whose junction is
//! already explained by the pangenome. A line is dropped if it lacks an `eseq`
//! tag, or if some ropebwt3 alignment's query interval spans (contains) the
//! junction interval on an eseq window (extended by `-l` on each side) — meaning
//! the breakend is "protected" and not novel. An insertion has two windows
//! (`eseq`, `eseq2`); it is dropped if either junction is spanned.
//!
//! The PAF's query names are `getfa`'s names (`readName_idx_leftFlank_
//! rightFlank_code`) in the same order as the eseq windows, and a subset of
//! them, so both files are streamed with a one-line lookahead over the PAF. A
//! PAF line left unconsumed at the end means the PAF was not made from this
//! input (in this order) and is an error, since later breakends would otherwise
//! pass as novel.
//!
//! With `-s STR` (`source: Some`), the mode flips: every input line is printed,
//! and the survivors (kept/novel breakends) get their `source=` INFO tag
//! rewritten to `STR`; dropped lines are printed verbatim with their original
//! source. This tags novel calls distinctly for a downstream `merge`.
//!
//! `-Q` (default 20) additionally drops any line whose `equal` quality (the
//! lower of `equal`/`equal2`) is below the threshold, in both modes (lines
//! without an `equal` tag are kept). A record without a known `ctype` is an
//! error.

use std::io::{self, BufRead, Write};

use crate::getfa::parse_windows;
use crate::io::{line_error, open_reader};

/// Filter `extract_out` using the ropebwt3 PAF `rb3_paf`; `l` is the margin.
///
/// With `source = Some(str)`, print all input lines and rewrite the `source=`
/// tag to `str` on survivors; with `None`, print only survivors, verbatim.
/// Lines whose `equal` quality is below `min_equal` are dropped in either mode.
pub fn run(
    extract_out: &str,
    rb3_paf: &str,
    l: i64,
    source: Option<&str>,
    min_equal: i64,
) -> io::Result<()> {
    let clip = open_reader(extract_out)?;
    let mut paf = open_reader(rb3_paf)?.lines();
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    // One-line lookahead into the PAF stream.
    let mut ahead: Option<String> = None;

    for (lineno, line) in clip.lines().enumerate() {
        let line = line?;
        let junctions =
            protected_intervals(&line, l).map_err(|e| line_error(extract_out, lineno + 1, e))?;
        let Some((junctions, equal)) = junctions.filter(|(j, _)| !j.is_empty()) else {
            // No eseq tag -> dropped. In print-all mode it is still emitted
            // verbatim; it has no PAF entry, so the stream stays in sync.
            if source.is_some() {
                writeln!(out, "{line}")?;
            }
            continue;
        };

        // Consume the PAF lines belonging to each eseq window (same order,
        // subset), and test whether any spans that window's [lo, hi]. Every
        // window is consumed, even once one is protected, to stay in sync.
        let mut protected = false;
        for (name, lo, hi) in &junctions {
            loop {
                if ahead.is_none() {
                    ahead = next_paf(&mut paf)?;
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
                    && qs <= *lo
                    && qe >= *hi
                {
                    protected = true;
                }
            }
        }

        // -Q: drop low-quality-eseq breakends (in either mode). Done AFTER the
        // PAF consumption above so the one-line lookahead stays in sync no matter
        // how the PAF was built.
        if equal.is_some_and(|q| q < min_equal) {
            continue;
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
    // Any PAF line left over matched no window in order: the PAF does not come
    // from `getfa` on this input, and the breakends after the stall were all
    // judged novel.
    if ahead.is_none() {
        ahead = next_paf(&mut paf)?;
    }
    if let Some(p) = ahead {
        let qn = p.split('\t').next().unwrap_or("");
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{rb3_paf}: PAF query `{qn}` matches no eseq window of {extract_out} in order; \
                 the PAF must come from `getfa` on the same input"
            ),
        ));
    }
    out.flush()?;
    Ok(())
}

/// The next non-empty PAF line, if any.
fn next_paf(paf: &mut impl Iterator<Item = io::Result<String>>) -> io::Result<Option<String>> {
    for p in paf {
        let p = p?;
        if !p.is_empty() {
            return Ok(Some(p));
        }
    }
    Ok(None)
}

/// Rewrite the `source=` value in a record's INFO (column 9) to `src`. If no
/// `source=` tag is present, prepend one (extract's "source first" convention).
/// Only this one tag's value changes; tag order is not significant.
fn relabel_source(line: &str, src: &str) -> String {
    let mut f: Vec<&str> = line.split('\t').collect();
    if f.len() < 9 {
        return line.to_string();
    }
    let mut found = false;
    let mut tags: Vec<String> = f[8]
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
    let info = tags.join(";");
    f[8] = &info;
    f.join("\t")
}

/// A protected junction interval `(fastaName, lo, hi)` on one eseq window.
type Junction = (String, i64, i64);
/// A line's junction intervals and its `equal` quality.
type Junctions = (Vec<Junction>, Option<i64>);

/// For an `extract` line, return the protected junction interval of each eseq
/// window (named as by `getfa`; empty if none), and the
/// line's `equal` quality (lowest over its windows). `Ok(None)` for a line with
/// fewer than 9 columns, treated like one without eseq; an error for a record
/// without a known `ctype`.
///
/// With `elen=e0,e1,e2`, the junction on the eseq is `[e0, e0+|e1|]`, extended
/// by `l` on each side and clamped to the eseq.
fn protected_intervals(line: &str, l: i64) -> Result<Option<Junctions>, String> {
    let Some(w) = parse_windows(line)? else {
        return Ok(None);
    };
    let junctions = w
        .windows
        .into_iter()
        .map(|win| {
            let [e0, e1, e2] = win.elen;
            let eseq_len = e0 + e1.abs() + e2;
            let lo = (e0 - l).max(0);
            let hi = (e0 + e1.abs() + l).min(eseq_len);
            (win.name, lo, hi)
        })
        .collect();
    Ok(Some((junctions, w.equal)))
}

/// Parse a PAF line's 0-based query interval `(qs, qe)` (columns 3 and 4).
fn paf_query_interval(line: &str) -> Option<(i64, i64)> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 4 {
        return None;
    }
    Some((f[2].parse().ok()?, f[3].parse().ok()?))
}
