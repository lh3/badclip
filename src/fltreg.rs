//! The `fltreg` subcommand: drop `extract`/`merge` breakends that fall in BED
//! regions.
//!
//! Both `extract` and `merge` output share the 9-column layout, so one filter
//! serves both. A line is dropped if **either** breakend endpoint — `(ctg1,pos1)`
//! (cols 0,1) or, when present, `(ctg2,pos2)` (cols 3,4) — lands inside a BED
//! region; clips (`ctg2="."`) only test the first endpoint. Survivors are printed
//! verbatim. Positions are raw 0-based offsets and BED intervals are half-open, so
//! an offset `p` is inside `[start,end)` iff `start <= p < end`.
//!
//! Regions are held per contig as start-sorted, merged (disjoint) intervals, so a
//! single binary search (`partition_point`) decides containment — no new
//! dependency.

use std::collections::HashMap;
use std::io::{self, BufRead, Write};

use crate::io::open_reader;

/// BED regions grouped by contig, each contig's intervals sorted by start and
/// merged so they are disjoint (see [`Regions::load`]).
struct Regions(HashMap<String, Vec<(i64, i64)>>);

impl Regions {
    /// Load a BED file (gzip ok; "-" for stdin) into per-contig merged intervals.
    /// Blank lines and `#`/`track`/`browser` header lines are skipped, as are
    /// lines whose first three tab columns are not `chrom start end` with integer
    /// coordinates.
    fn load(path: &str) -> io::Result<Regions> {
        let mut map: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for line in open_reader(path)?.lines() {
            let line = line?;
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }
            let mut f = line.split('\t');
            let (Some(chrom), Some(s), Some(e)) = (f.next(), f.next(), f.next()) else {
                continue;
            };
            if let (Ok(s), Ok(e)) = (s.parse::<i64>(), e.parse::<i64>()) {
                map.entry(chrom.to_string()).or_default().push((s, e));
            }
        }
        // Sort by start and merge overlapping/adjacent intervals so each contig's
        // list is disjoint and start-sorted — the precondition for the binary
        // search in `contains` to be correct even on an overlapping BED.
        for v in map.values_mut() {
            v.sort_unstable();
            let mut merged: Vec<(i64, i64)> = Vec::with_capacity(v.len());
            for &(s, e) in v.iter() {
                match merged.last_mut() {
                    Some(last) if s <= last.1 => last.1 = last.1.max(e),
                    _ => merged.push((s, e)),
                }
            }
            *v = merged;
        }
        Ok(Regions(map))
    }

    /// Whether offset `pos` falls in some region on `ctg` (half-open `[s,e)`).
    fn contains(&self, ctg: &str, pos: i64) -> bool {
        let Some(v) = self.0.get(ctg) else {
            return false;
        };
        // Rightmost interval with start <= pos; disjoint+sorted, so it is the only
        // one that could contain pos.
        let i = v.partition_point(|iv| iv.0 <= pos);
        i > 0 && pos < v[i - 1].1
    }
}

/// Filter `input` (`extract`/`merge` output; gzip ok, "-" for stdin) against the
/// BED file `bed`, printing lines whose breakends both avoid every region.
pub fn run(input: &str, bed: &str) -> io::Result<()> {
    let regions = Regions::load(bed)?;
    let reader = open_reader(input)?;
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    for line in reader.lines() {
        let line = line?;
        let mut f = line.split('\t');
        // Endpoint 1 = cols 0,1; endpoint 2 = cols 3,4 (absent for a clip).
        let ctg1 = f.next();
        let pos1 = f.next().and_then(|s| s.parse::<i64>().ok());
        let ctg2 = f.nth(1); // skip col 2 (ori)
        let pos2 = f.next().and_then(|s| s.parse::<i64>().ok());

        let hit1 = matches!((ctg1, pos1), (Some(c), Some(p)) if regions.contains(c, p));
        let hit2 = matches!((ctg2, pos2), (Some(c), Some(p)) if regions.contains(c, p));
        if !hit1 && !hit2 {
            writeln!(out, "{line}")?;
        }
    }
    Ok(())
}
