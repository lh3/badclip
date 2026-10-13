//! The `getfa` subcommand (alias `geteseq`): turn `extract` output into a FASTA
//! of the extracted breakend sequences.
//!
//! Reads `extract`'s TAB-delimited records and, for each eseq window a record
//! carries, writes a FASTA record whose name is
//! `readName_idx_leftFlank_rightFlank_code` (from the `idx` and `elen` tags and
//! the breakend type) and whose sequence is the window's bases. The code is `T`
//! (trim), `S` (split), `D` (del), or, for an insertion's two windows, `L`
//! (`eseq`, its first breakpoint along the read) and `R` (`eseq2`, the second).
//! Records without an eseq are skipped, as are records whose `equal` quality
//! (the lower of `equal`/`equal2`) is below `-Q` (records lacking `equal` are
//! kept). The parsing here is shared with `fltbnd`, so names always agree.

use std::io::{self, BufRead, Write};

use crate::io::open_reader;

/// One eseq window of an `extract` record.
pub(crate) struct Window<'a> {
    /// FASTA name: `readName_idx_leftFlank_rightFlank_code`.
    pub name: String,
    /// `elen=leftFlank,qdist,rightFlank` of the window.
    pub elen: [i64; 3],
    pub eseq: &'a str,
}

/// The eseq windows of an `extract` record, in FASTA order, and its quality.
pub(crate) struct Windows<'a> {
    /// Windows that carry an eseq (two for an insertion, else at most one).
    pub windows: Vec<Window<'a>>,
    /// The lowest `equal`/`equal2` value, if any.
    pub equal: Option<i64>,
}

/// Parse the eseq windows of one `extract` record. `None` for a malformed line
/// or an unknown `ctype`. A record without `ctype` (pre-0.2 output) is a trim if
/// it has no mate (`ctg2` is `.`), else a split.
pub(crate) fn parse_windows(line: &str) -> Option<Windows<'_>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }
    let qname = fields[5];

    let mut ctype = None;
    let mut idx = None;
    // Index 0: elen/eseq; index 1: elen2/eseq2 (an insertion's second window).
    let mut elen = [None, None];
    let mut eseq = [None, None];
    let mut equal: Option<i64> = None;
    for kv in fields[8].split(';') {
        let Some((k, v)) = kv.split_once('=') else {
            continue;
        };
        match k {
            "ctype" => ctype = Some(v),
            "idx" => idx = Some(v),
            "elen" => elen[0] = Some(v),
            "elen2" => elen[1] = Some(v),
            "eseq" => eseq[0] = Some(v),
            "eseq2" => eseq[1] = Some(v),
            "equal" | "equal2" => {
                if let Ok(q) = v.parse::<i64>() {
                    equal = Some(equal.map_or(q, |e| e.min(q)));
                }
            }
            _ => {}
        }
    }
    let codes = match ctype {
        Some("trim") => ["T", ""],
        Some("split") => ["S", ""],
        Some("del") => ["D", ""],
        Some("ins") => ["L", "R"],
        Some(_) => return None,
        None if fields[3] == "." => ["T", ""],
        None => ["S", ""],
    };

    let mut windows = Vec::new();
    for i in 0..2 {
        let (Some(elen), Some(eseq)) = (elen[i], eseq[i]) else {
            continue;
        };
        if codes[i].is_empty() {
            continue;
        }
        let e: Vec<i64> = elen.split(',').filter_map(|x| x.parse().ok()).collect();
        if e.len() < 3 {
            continue;
        }
        let elen = [e[0], e[1], e[2]];
        windows.push(Window {
            name: format!("{qname}_{}_{}_{}_{}", idx?, elen[0], elen[2], codes[i]),
            elen,
            eseq,
        });
    }
    Some(Windows { windows, equal })
}

/// Read `extract` output from `input` (`"-"` = stdin) and emit FASTA. Records
/// whose `equal` quality is below `min_equal` are dropped.
pub fn run(input: &str, min_equal: i64) -> io::Result<()> {
    let reader = open_reader(input)?;
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    for line in reader.lines() {
        let line = line?;
        let Some(w) = parse_windows(&line) else {
            continue;
        };
        // Drop records with a low eseq quality; keep those lacking an `equal` tag.
        if w.equal.is_some_and(|q| q < min_equal) {
            continue;
        }
        for win in &w.windows {
            writeln!(out, ">{}", win.name)?;
            writeln!(out, "{}", win.eseq)?;
        }
    }
    Ok(())
}
