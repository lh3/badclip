//! PAF parsing.
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
}

/// Parse a single PAF line into a [`Hit`].
///
/// Returns `None` when the line is not a usable primary/supplementary
/// alignment: too few columns, a non-`+/-` strand field, unparseable
/// coordinates, or an alignment-type tag (`tp:A:`) that is not `P`.
pub fn parse_paf_line(line: &str) -> Option<Hit> {
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

    // Keep only primary/supplementary alignments (tp:A:P). minimap2 marks both
    // primary and supplementary chimeric segments with tp:A:P; secondary
    // alignments carry tp:A:S. Lines without a tp tag are dropped as well,
    // matching the reference behaviour.
    let mut is_primary = false;
    for tag in &t[12..] {
        if let Some(rest) = tag.strip_prefix("tp:A:") {
            is_primary = rest == "P";
            break;
        }
    }
    if !is_primary {
        return None;
    }

    Some(Hit {
        qname: t[0].to_string(),
        qlen: t[1].parse().ok()?,
        qs: t[2].parse().ok()?,
        qe: t[3].parse().ok()?,
        strand,
        ctg: t[5].to_string(),
        ts: t[7].parse().ok()?,
        te: t[8].parse().ok()?,
        alen: t[10].parse().ok()?,
        mapq: t[11].parse().ok()?,
    })
}
