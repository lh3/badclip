//! The `fltcnt` subcommand: keep `merge` calls with enough raw/filtered support.
//!
//! A `merge` line's `count=` INFO tag lists, per `source=`, the forward,reverse
//! read counts (`count=S1:f,r|S2:f,r`). The expected input merges two kinds of
//! source: the raw breakends (named `*.raw`) and the `flteseq -s` survivors
//! (named `*.flt`). `--src-raw`/`--src-flt` replace the suffix match with an exact
//! source name. Counts of every matching source are summed. A line with neither
//! a raw nor a filtered source is an error (the naming convention doesn't fit,
//! so the sources must be given explicitly).
//!
//! A line is printed verbatim if the filtered reads total `>= -c`, each strand
//! of the filtered reads has `>= -s`, and the raw reads total `>= -r`. Unless
//! `-l` is negative, clips (`ctg2 = "."`) and same-contig calls with
//! `pos2 - pos1 < -l` are also dropped (so the default `-l 0` drops clips only);
//! inter-contig joins have no distance and pass.

use std::io::{self, BufRead, Write};

use crate::io::open_reader;

pub struct FltcntOpts {
    pub input: String,
    pub min_flt: i64,
    pub min_flt_strand: i64,
    pub min_raw: i64,
    pub min_dist: i64,
    pub src_raw: Option<String>,
    pub src_flt: Option<String>,
}

/// Summed `(fwd, rev)` counts of the raw and filtered sources of one line, or
/// `None` if the line has no `count=` tag at all.
struct Counts {
    raw: Option<(i64, i64)>,
    flt: Option<(i64, i64)>,
}

/// Whether `src` names a source of the given kind: exact match with the explicit
/// name if one was given, else the `suffix` convention.
fn is_src(src: &str, name: Option<&str>, suffix: &str) -> bool {
    match name {
        Some(n) => src == n,
        None => src.ends_with(suffix),
    }
}

fn add(acc: &mut Option<(i64, i64)>, f: i64, r: i64) {
    let a = acc.get_or_insert((0, 0));
    a.0 += f;
    a.1 += r;
}

/// Parse the `count=` tag of an INFO column into raw/filtered counts.
fn parse_counts(info: &str, o: &FltcntOpts) -> Option<Counts> {
    let val = info.split(';').find_map(|t| t.strip_prefix("count="))?;
    let mut c = Counts { raw: None, flt: None };
    for ent in val.split('|') {
        // An entry without a `src:` prefix (e.g. `merge -m` output) matches
        // neither kind.
        let Some((src, fr)) = ent.rsplit_once(':') else {
            continue;
        };
        let Some((f, r)) = fr.split_once(',') else {
            continue;
        };
        let (Ok(f), Ok(r)) = (f.parse::<i64>(), r.parse::<i64>()) else {
            continue;
        };
        if is_src(src, o.src_raw.as_deref(), ".raw") {
            add(&mut c.raw, f, r);
        }
        if is_src(src, o.src_flt.as_deref(), ".flt") {
            add(&mut c.flt, f, r);
        }
    }
    Some(c)
}

pub fn run(o: &FltcntOpts) -> io::Result<()> {
    let reader = open_reader(&o.input)?;
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let counts = f.get(8).and_then(|info| parse_counts(info, o));
        let Some(Counts { raw, flt }) = counts.filter(|c| c.raw.is_some() || c.flt.is_some())
        else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "line {}: no raw (`{}`) or filtered (`{}`) source in `count=`; \
                     specify the sources with --src-raw and --src-flt",
                    lineno + 1,
                    o.src_raw.as_deref().unwrap_or("*.raw"),
                    o.src_flt.as_deref().unwrap_or("*.flt"),
                ),
            ));
        };
        let (rf, rr) = raw.unwrap_or((0, 0));
        let (ff, fr) = flt.unwrap_or((0, 0));
        if ff + fr < o.min_flt || ff.min(fr) < o.min_flt_strand || rf + rr < o.min_raw {
            continue;
        }
        if o.min_dist >= 0 {
            // Clip: no second endpoint.
            if f.get(3).is_none_or(|&c| c == ".") {
                continue;
            }
            if o.min_dist > 0 && f.first() == f.get(3) {
                let pos = |i: usize| f.get(i).and_then(|s| s.parse::<i64>().ok());
                if let (Some(p1), Some(p2)) = (pos(1), pos(4)) {
                    if p2 - p1 < o.min_dist {
                        continue;
                    }
                }
            }
        }
        writeln!(out, "{line}")?;
    }
    Ok(())
}
