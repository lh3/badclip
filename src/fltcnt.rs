//! The `fltcnt` subcommand: keep `merge` calls with enough read support from
//! selected sources.
//!
//! A `merge` line's `count=` INFO tag lists, per `source=`, the forward,reverse
//! read counts (`count=S1:f,r|S2:f,r`). Sources fall into two groups: `--src`
//! (comma-separated names; default: every source ending with `.flt`, i.e. the
//! `flteseq -s` survivors) and `--rest` (comma-separated names; default: every
//! source not in `--src`). Each group's counts are summed over its sources. A
//! line with no source in either group (e.g. `merge -m` output, whose `count=`
//! has no source labels) is an error.
//!
//! A line is printed verbatim if the `--src` reads total `>= -c`, each strand of
//! the `--src` reads has `>= -s`, and the `--rest` reads total `>= -r`. Unless
//! `-l` is negative, clips (`ctg2 = "."`) and same-contig calls with
//! `pos2 - pos1 < -l` are also dropped (so the default `-l 0` drops clips only);
//! inter-contig joins have no distance and pass.

use std::io::{self, BufRead, Write};

use crate::io::open_reader;

pub struct FltcntOpts {
    pub input: String,
    pub min_src: i64,
    pub min_src_strand: i64,
    pub min_rest: i64,
    pub min_dist: i64,
    /// `--src` names; empty = sources ending with `.flt`.
    pub src: Vec<String>,
    /// `--rest` names; empty = sources not in `--src`.
    pub rest: Vec<String>,
}

impl FltcntOpts {
    fn is_src(&self, s: &str) -> bool {
        if self.src.is_empty() {
            s.ends_with(".flt")
        } else {
            self.src.iter().any(|n| n == s)
        }
    }

    fn is_rest(&self, s: &str) -> bool {
        if self.rest.is_empty() {
            !self.is_src(s)
        } else {
            self.rest.iter().any(|n| n == s)
        }
    }
}

/// Summed `(fwd, rev)` counts of the `--src` and `--rest` sources of one line;
/// `None` for a group with no source on the line.
struct Counts {
    src: Option<(i64, i64)>,
    rest: Option<(i64, i64)>,
}

fn add(acc: &mut Option<(i64, i64)>, f: i64, r: i64) {
    let a = acc.get_or_insert((0, 0));
    a.0 += f;
    a.1 += r;
}

/// Parse the `count=` tag of an INFO column into `--src`/`--rest` counts, or
/// `None` if the line has no `count=` tag at all.
fn parse_counts(info: &str, o: &FltcntOpts) -> Option<Counts> {
    let val = info.split(';').find_map(|t| t.strip_prefix("count="))?;
    let mut c = Counts { src: None, rest: None };
    for ent in val.split('|') {
        // An entry without a `src:` prefix (e.g. `merge -m` output) matches
        // neither group.
        let Some((src, fr)) = ent.rsplit_once(':') else {
            continue;
        };
        let Some((f, r)) = fr.split_once(',') else {
            continue;
        };
        let (Ok(f), Ok(r)) = (f.parse::<i64>(), r.parse::<i64>()) else {
            continue;
        };
        if o.is_src(src) {
            add(&mut c.src, f, r);
        }
        if o.is_rest(src) {
            add(&mut c.rest, f, r);
        }
    }
    Some(c)
}

pub fn run(o: &FltcntOpts) -> io::Result<()> {
    let reader = open_reader(&o.input)?;
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    // Whether any line had a `--src` source, to warn when none does (e.g. the
    // default `*.flt` naming doesn't fit the input).
    let mut seen_src = false;
    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let counts = f.get(8).and_then(|info| parse_counts(info, o));
        let Some(Counts { src, rest }) = counts.filter(|c| c.src.is_some() || c.rest.is_some())
        else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "line {}: no --src or --rest source in `count=`; \
                     specify the sources with --src and --rest",
                    lineno + 1,
                ),
            ));
        };
        seen_src |= src.is_some();
        let (sf, sr) = src.unwrap_or((0, 0));
        let (rf, rr) = rest.unwrap_or((0, 0));
        if sf + sr < o.min_src || sf.min(sr) < o.min_src_strand || rf + rr < o.min_rest {
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
    if !seen_src {
        let src = if o.src.is_empty() { "*.flt".to_string() } else { o.src.join(",") };
        eprintln!("badclip: warning: no line has a --src source ({src}); specify it with --src");
    }
    Ok(())
}
