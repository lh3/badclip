//! PAF parsing, plus the CIGAR helpers shared with the alignment-file path.
//!
//! A [`Hit`] is one primary/supplementary alignment of a read to a contig,
//! reduced to just the fields the breakend logic needs. Coordinates use the
//! project's offset convention: an interval `|st,en|` marks offsets *between*
//! bases, i.e. the usual 0-based half-open range.

/// Strand of an alignment.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Strand {
    Fwd,
    Rev,
}

impl Strand {
    pub fn as_char(self) -> char {
        match self {
            Strand::Fwd => '+',
            Strand::Rev => '-',
        }
    }
}

/// A single alignment (hit) of a read/query onto a target contig.
///
/// The query interval `|qs,qe|` maps to `ctg:|ts,te|` on `strand`.
#[derive(Clone, Debug)]
pub struct Hit {
    pub qname: String,
    pub qlen: i64,
    pub qs: i64,
    pub qe: i64,
    pub strand: Strand,
    pub ctg: String,
    pub ts: i64,
    pub te: i64,
    /// Alignment block length (PAF column 11).
    pub alen: i64,
    pub mapq: i64,
    /// Long indels in this hit's CIGAR (`-g`). Empty when the full CIGAR is not
    /// known, e.g. for a hit rebuilt from a (collapsed) `SA` tag.
    pub gaps: Vec<Gap>,
}

/// A long indel (`I`/`D` op of at least `-g` bases) inside one alignment's
/// CIGAR, in alignment (reference-forward) orientation.
#[derive(Clone, Debug)]
pub struct Gap {
    /// An insertion (`I`); otherwise a deletion (`D`).
    pub ins: bool,
    /// Reference offset of the CIGAR placement: a deletion removes
    /// `|tpos,tpos+len|`; an insertion sits at offset `tpos`.
    pub tpos: i64,
    /// Leftmost and rightmost reference offsets at which an insertion can be
    /// placed equivalently (see [`find_gaps`]). Both equal `tpos` for a deletion
    /// or when the read sequence is unavailable.
    pub tpos_l: i64,
    pub tpos_r: i64,
    /// Indel length.
    pub len: i64,
    /// Query offset of the event in alignment orientation, counting any leading
    /// clip (for an insertion, the offset of its first inserted base).
    pub qoff: i64,
}

/// Parse a CIGAR string into `(op_char, len)` pairs.
pub fn parse_cigar(s: &str) -> Vec<(char, i64)> {
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

fn is_match_op(op: char) -> bool {
    matches!(op, 'M' | '=' | 'X')
}

/// Collect the indels of at least `min_gap` bases (`min_gap <= 0` disables
/// this) from a CIGAR. The walk starts at query offset `q0` and reference
/// offset `t0`; clip ops (`S`/`H`) consume query, so a full SAM CIGAR starts at
/// `q0 = 0` while a PAF `cg:Z` (aligned part only) starts at the aligned start.
///
/// `seq`, when given, is the stored read sequence (alignment orientation) with
/// the query offset of its first base (the leading hard-clip length). It is used
/// to normalize each insertion: with the inserted bases at `seq[a..b]`, the
/// insertion can move right while `seq[a+i] == seq[b+i]` and left while
/// `seq[a-1-i] == seq[b-1-i]`. A shift stays inside the run of `M`/`=`/`X` ops
/// next to the insertion on that side (so every shifted base is one reference
/// base) and stops at an `N` or the end of `seq`.
pub fn find_gaps(
    ops: &[(char, i64)],
    q0: i64,
    t0: i64,
    min_gap: i64,
    seq: Option<(&[u8], i64)>,
) -> Vec<Gap> {
    let mut gaps = Vec::new();
    if min_gap <= 0 {
        return gaps;
    }
    let (mut q, mut t) = (q0, t0);
    for (i, &(op, len)) in ops.iter().enumerate() {
        match op {
            'M' | '=' | 'X' => {
                q += len;
                t += len;
            }
            'I' => {
                if len >= min_gap {
                    let (kl, kr) = match seq {
                        Some((s, s0)) => {
                            // Lengths of the match runs flanking this op.
                            let run = |it: &mut dyn Iterator<Item = &(char, i64)>| -> i64 {
                                it.take_while(|(o, _)| is_match_op(*o)).map(|(_, l)| l).sum()
                            };
                            let max_l = run(&mut ops[..i].iter().rev());
                            let max_r = run(&mut ops[i + 1..].iter());
                            shift_ins(s, q - s0, len, max_l, max_r)
                        }
                        None => (0, 0),
                    };
                    gaps.push(Gap {
                        ins: true,
                        tpos: t,
                        tpos_l: t - kl,
                        tpos_r: t + kr,
                        len,
                        qoff: q,
                    });
                }
                q += len;
            }
            'D' => {
                if len >= min_gap {
                    gaps.push(Gap {
                        ins: false,
                        tpos: t,
                        tpos_l: t,
                        tpos_r: t,
                        len,
                        qoff: q,
                    });
                }
                t += len;
            }
            'N' => t += len,
            'S' | 'H' => q += len,
            _ => {} // P (pad)
        }
    }
    gaps
}

/// How far `(left, right)` an insertion of `len` bases at `seq[a..a+len]` can
/// be shifted, at most `max_l`/`max_r` bases (see [`find_gaps`]).
fn shift_ins(seq: &[u8], a: i64, len: i64, max_l: i64, max_r: i64) -> (i64, i64) {
    let n = seq.len() as i64;
    let b = a + len;
    if a < 0 || b > n {
        return (0, 0); // seq does not cover the insertion
    }
    let same = |x: u8, y: u8| {
        let (x, y) = (x.to_ascii_uppercase(), y.to_ascii_uppercase());
        x == y && x != b'N'
    };
    let mut kr = 0;
    while kr < max_r && b + kr < n && same(seq[(a + kr) as usize], seq[(b + kr) as usize]) {
        kr += 1;
    }
    let mut kl = 0;
    while kl < max_l && a - 1 - kl >= 0 && same(seq[(a - 1 - kl) as usize], seq[(b - 1 - kl) as usize]) {
        kl += 1;
    }
    (kl, kr)
}

/// If `line` is an unmapped-read PAF line (target name `*`, as minimap2
/// writes with `--paf-no-hit`), return its query length.
pub fn parse_paf_unmapped(line: &str) -> Option<i64> {
    let mut t = line.split('\t');
    let _qname = t.next()?;
    let qlen = t.next()?.parse().ok()?;
    (t.nth(3)? == "*").then_some(qlen)
}

/// Parse a single PAF line into a [`Hit`]. Indels of at least `min_gap` bases
/// are collected from the `cg:Z:` tag when present (no read sequence, so
/// insertions are not normalized).
///
/// Returns `None` when the line is not a usable alignment: too few columns, a
/// non-`+/-` strand field, unparseable coordinates, or an alignment-type tag
/// (`tp:A:`) that is neither `P` (primary/supplementary) nor `I` (inversion).
pub fn parse_paf_line(line: &str, min_gap: i64) -> Option<Hit> {
    let t: Vec<&str> = line.split('\t').collect();
    // PAF has at least 12 mandatory columns.
    if t.len() < 12 {
        return None;
    }

    let strand = match t[4] {
        "+" => Strand::Fwd,
        "-" => Strand::Rev,
        _ => return None,
    };

    // Keep primary/supplementary (tp:A:P) and inversion (tp:A:I) alignments.
    // minimap2 marks primary and supplementary chimeric segments with tp:A:P and
    // the inverted middle segment of an inversion with tp:A:I; both are real
    // hits for breakend calling (and BAM exposes the inversion via the SA tag).
    // Secondary (tp:A:S) and tp-less lines are dropped.
    let mut keep = false;
    let mut cigar = None;
    for tag in &t[12..] {
        if let Some(rest) = tag.strip_prefix("tp:A:") {
            keep = rest == "P" || rest == "I";
        } else if let Some(cg) = tag.strip_prefix("cg:Z:") {
            cigar = Some(cg);
        }
    }
    if !keep {
        return None;
    }

    let qlen: i64 = t[1].parse().ok()?;
    let qs: i64 = t[2].parse().ok()?;
    let qe: i64 = t[3].parse().ok()?;
    let ts: i64 = t[7].parse().ok()?;
    // `cg:Z` covers the aligned part only; in alignment orientation it starts
    // at query offset `qs` (forward) or `qlen - qe` (reverse).
    let gaps = match cigar {
        Some(cg) => {
            let q0 = if strand == Strand::Fwd { qs } else { qlen - qe };
            find_gaps(&parse_cigar(cg), q0, ts, min_gap, None)
        }
        None => Vec::new(),
    };
    Some(Hit {
        qname: t[0].to_string(),
        qlen,
        qs,
        qe,
        strand,
        ctg: t[5].to_string(),
        ts,
        te: t[8].parse().ok()?,
        alen: t[10].parse().ok()?,
        mapq: t[11].parse().ok()?,
        gaps,
    })
}

#[cfg(test)]
mod tests {
    use super::{find_gaps, parse_cigar};

    /// `(ins, tpos, tpos_l, tpos_r, len, qoff)` of each gap.
    fn gaps(cigar: &str, min_gap: i64, seq: Option<(&[u8], i64)>) -> Vec<(bool, i64, i64, i64, i64, i64)> {
        find_gaps(&parse_cigar(cigar), 0, 100, min_gap, seq)
            .into_iter()
            .map(|g| (g.ins, g.tpos, g.tpos_l, g.tpos_r, g.len, g.qoff))
            .collect()
    }

    #[test]
    fn finds_long_indels() {
        // 2S + 10M, then a 5D and a 4I; 1-bp indels stay below -g 3.
        let g = gaps("2S10M5D3M1I3M4I6M1D2M", 3, None);
        assert_eq!(g, vec![(false, 110, 110, 110, 5, 12), (true, 121, 121, 121, 4, 19)]);
        assert!(gaps("2S10M5D3M", 0, None).is_empty(), "-g 0 disables");
    }

    #[test]
    fn insertion_shifts() {
        // seq: AAACG | TTCG (inserted) | TTCGA ; the insertion repeats the next
        // 4 read bases (TSD-like), then 'A' vs 'T' stops the right shift. On the
        // left, seq[a-1]='G' == seq[b-1]='G', 'C'=='C', 'A' vs 'T' stops -> 2.
        let seq = b"AAACGTTCGTTCGA";
        let g = gaps("5M4I5M", 4, Some((seq, 0)));
        assert_eq!(g, vec![(true, 105, 103, 109, 4, 5)]);
        // A leading hard clip shifts the SEQ index: same read, 3H prepended.
        let g = gaps("3H5M4I5M", 4, Some((seq, 3)));
        assert_eq!(g, vec![(true, 105, 103, 109, 4, 8)]);
    }

    #[test]
    fn insertion_shift_bounds() {
        let seq = b"AAACGTTCGTTCGA";
        // The right match run is only 2 bases (then a deletion) -> kr <= 2; the
        // left run 1 base (after a 1-bp insertion) -> kl <= 1.
        let g = gaps("3M1I1M4I2M1D3M", 4, Some((seq, 0)));
        assert_eq!(g, vec![(true, 104, 103, 106, 4, 5)]);
        // N never matches.
        let seq = b"AAACGNTCGNTCGA";
        let g = gaps("5M4I5M", 4, Some((seq, 0)));
        assert_eq!(g, vec![(true, 105, 103, 105, 4, 5)]);
        // A tandem repeat: the shift may exceed the insertion length.
        let seq = b"GCACACACACT";
        let g = gaps("3M2I6M", 2, Some((seq, 0)));
        assert_eq!(g, vec![(true, 103, 101, 108, 2, 3)]);
    }
}
