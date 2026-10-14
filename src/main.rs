//! badclip — extract breakends and structural-variant signals from long-read
//! alignments.

mod aln;
mod extract;
mod fltbnd;
mod fltcnt;
mod fltreg;
mod getfa;
#[allow(dead_code)]
mod iitree;
mod io;
mod merge;
mod paf;
mod record;

use std::io::Write;
use std::process::ExitCode;

use clap::{CommandFactory, Parser, Subcommand};

use extract::ExtractOpts;
use fltcnt::FltcntOpts;
use merge::MergeOpts;

#[derive(Parser)]
#[command(name = "badclip", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Extract breakends from read alignments in SAM/BAM/CRAM (or PAF) format.
    Extract {
        /// Input alignments: SAM/BAM/CRAM (autodetected) by default, or PAF with
        /// --paf. Use "-" for stdin.
        input: Option<String>,

        /// Read PAF (optionally gzip'd) instead of SAM/BAM/CRAM.
        #[arg(long)]
        paf: bool,

        /// Reference FASTA (faidx-indexed), required for CRAM input.
        #[arg(short = 'r', long = "reference")]
        reference: Option<String>,

        /// File listing ALT contig names; reads/SA hits on them are excluded.
        #[arg(long = "alt", value_name = "FILE")]
        alt: Option<String>,

        /// Dataset name, emitted as a `source=` INFO tag on every record.
        #[arg(short = 's', long = "source", default_value_t = String::from("foo"))]
        source: String,

        /// Report a trim breakend only when the clip is longer than this (k/m/g
        /// suffix ok).
        #[arg(short = 'c', long = "min-clip", default_value_t = 50, value_parser = parse_num)]
        min_clip: i64,

        /// Also report trim breakends (read ends clipped away from every
        /// alignment).
        #[arg(short = 'T', long = "trim")]
        trim: bool,

        /// Minimum CIGAR indel length to report a del/ins breakend (0 disables;
        /// k/m/g suffix ok).
        #[arg(short = 'g', long = "min-gap", default_value_t = 1000, value_parser = parse_num)]
        min_gap: i64,

        /// Drop hits whose alignment block length (PAF col 11) is below this
        /// (k/m/g suffix ok).
        #[arg(short = 'a', long = "min-aln-len", default_value_t = 0, value_parser = parse_num)]
        min_aln_len: i64,

        /// Flanking read sequence extracted on each side of a breakend (not for
        /// PAF; k/m/g suffix ok).
        #[arg(short = 'f', long = "flank", default_value_t = 250, value_parser = parse_num)]
        flank: i64,

        /// Maximum extracted sequence length; longer windows omit eseq (not for
        /// PAF; k/m/g suffix ok).
        #[arg(short = 'e', long = "max-eseq", default_value_t = 1000, value_parser = parse_num)]
        max_eseq: i64,

        /// Drop breakends whose eseq quality (`equal`) is below this (0 = keep all).
        #[arg(short = 'Q', long = "min-equal", default_value_t = 0)]
        min_equal: i64,

        /// Drop output lines whose col-7 mapq is below this (post-filter).
        #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
        min_mapq: i64,
    },

    /// Convert `extract` output into a FASTA of extracted breakend sequences.
    #[command(alias = "geteseq")]
    Getfa {
        /// Drop records whose eseq quality (`equal` tag) is below this.
        #[arg(short = 'Q', long = "min-equal", default_value_t = 20)]
        min_equal: i64,

        /// `extract` output (gzip ok; "-" for stdin).
        input: Option<String>,
    },

    /// Filter `extract` breakends whose eseq is spanned by a pangenome alignment.
    #[command(alias = "flteseq")]
    Fltbnd {
        /// Margin extending the protected junction interval on each side (k/m/g
        /// suffix ok).
        #[arg(short = 'l', long = "margin", default_value_t = 50, value_parser = parse_num)]
        margin: i64,

        /// Print all input lines (except those dropped by -Q), rewriting
        /// survivors' `source=` tag to this.
        #[arg(short = 's', long = "source")]
        source: Option<String>,

        /// Drop lines whose eseq quality (`equal` tag) is below this.
        #[arg(short = 'Q', long = "min-equal", default_value_t = 20)]
        min_equal: i64,

        /// `extract` output (gzip ok; "-" for stdin).
        extract_out: Option<String>,

        /// ropebwt3 `sw` PAF on the `getfa` FASTA (gzip ok).
        rb3_paf: Option<String>,
    },

    /// Merge per-read `extract` breakends into consensus SV calls.
    Merge {
        /// `extract` output files, merged as if concatenated (gzip ok; "-" for
        /// stdin).
        inputs: Vec<String>,

        /// Minimum read count to emit a call.
        #[arg(short = 'c', long = "min-cnt", default_value_t = 2)]
        min_cnt: i64,

        /// Minimum read count on each strand.
        #[arg(short = 's', long = "min-cnt-strand", default_value_t = 0)]
        min_cnt_strand: i64,

        /// Clustering window size (bp; k/m/g suffix ok).
        #[arg(short = 'w', long = "win-size", default_value_t = 100, value_parser = parse_num)]
        win_size: i64,

        /// Cap on active clusters (flush trigger).
        #[arg(short = 'A', long = "max-allele", default_value_t = 100)]
        max_allele: i64,

        /// Maximum reads compared per cluster (deterministic cap).
        #[arg(short = 'M', long = "max-check", default_value_t = 500)]
        max_check: i64,

        /// Drop input breakends whose eseq quality (`equal` tag) is below this.
        #[arg(short = 'Q', long = "min-equal", default_value_t = 20)]
        min_equal: i64,

        /// Drop input breakends whose col-7 mapq (with -m: the higher of
        /// `avg_mapq`) is below this.
        #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
        min_mapq: i64,

        /// Drop input join breakends whose smaller per-hit mapq is below this
        /// (no effect on clips).
        #[arg(short = 'p', long = "min-mapq-min", default_value_t = 10)]
        min_mapq_min: i64,

        /// Treat input as combined `merge` output (sample-merge) instead of
        /// `extract` output.
        #[arg(short = 'm', long = "merge-input")]
        merge_input: bool,

        /// (-m only) drop input lines whose total count is below this.
        #[arg(short = 'C', long = "min-cnt-in", default_value_t = 3)]
        min_cnt_in: i64,

        /// (-m only) drop input lines whose per-strand count is below this.
        #[arg(short = 'S', long = "min-cnt-strand-in", default_value_t = 1)]
        min_cnt_strand_in: i64,
    },

    /// Filter `extract`/`merge` breakends that fall in BED regions.
    Fltreg {
        /// `extract`/`merge` output (gzip ok; "-" for stdin).
        input: Option<String>,

        /// BED file of regions to filter against (gzip ok).
        bed: Option<String>,

        /// Pad each BED interval to [st-l, en+l) before testing (k/m/g suffix
        /// ok).
        #[arg(short = 'l', long = "margin", default_value_t = 0, value_parser = parse_num)]
        margin: i64,
    },

    /// Keep `merge` calls with enough read support from selected sources.
    Fltcnt {
        /// `merge` output (gzip ok; "-" for stdin).
        input: Option<String>,

        /// Minimum reads from the --src sources.
        #[arg(short = 'c', long = "min-src", default_value_t = 3)]
        min_src: i64,

        /// Minimum reads from the --src sources on each strand.
        #[arg(short = 's', long = "min-src-strand", default_value_t = 0)]
        min_src_strand: i64,

        /// Maximum reads from the --rest sources (-1 = no limit).
        #[arg(
            short = 'r',
            long = "max-rest",
            default_value_t = 0,
            allow_negative_numbers = true
        )]
        max_rest: i64,

        /// Minimum pos2-pos1 distance (k/m/g suffix ok, e.g. 100k); any value
        /// >= 0 also drops clips, -1 disables the filter.
        #[arg(short = 'l', long = "min-dist", default_value_t = 0,
              value_parser = parse_num, allow_negative_numbers = true)]
        min_dist: i64,

        /// Apply -l to inversions (ori >< or <>) too [default: inversions pass
        /// regardless of distance].
        #[arg(short = 'I', long = "inv-dist")]
        inv_dist: bool,

        /// Maximum ins_nearby reads, summed over sources (-1 = no limit).
        #[arg(
            short = 'i',
            long = "max-ins-nearby",
            default_value_t = 0,
            allow_negative_numbers = true
        )]
        max_ins: i64,

        /// Comma-separated sources tested by -c/-s [default: sources ending
        /// with ".flt"].
        #[arg(long = "src", value_delimiter = ',')]
        src: Vec<String>,

        /// Comma-separated sources tested by -r [default: sources not in --src].
        #[arg(long = "rest", value_delimiter = ',')]
        rest: Vec<String>,
    },
}

/// Parse an integer with an optional case-insensitive `k`/`m`/`g` suffix
/// (×1e3/1e6/1e9), e.g. `100k`, `1m`, `1.5k`.
fn parse_num(s: &str) -> Result<i64, String> {
    let (num, mul) = match s.as_bytes().last().map(u8::to_ascii_lowercase) {
        Some(b'k') => (&s[..s.len() - 1], 1e3),
        Some(b'm') => (&s[..s.len() - 1], 1e6),
        Some(b'g') => (&s[..s.len() - 1], 1e9),
        _ => return s.parse::<i64>().map_err(|e| e.to_string()),
    };
    let x: f64 = num.parse().map_err(|_| format!("invalid number `{s}`"))?;
    Ok((x * mul).round() as i64)
}

/// Print a subcommand's help (as for `-h`) and return exit code 2. Used when a
/// subcommand is invoked with no input, so it shows usage instead of blocking on
/// stdin or erroring.
fn print_subcommand_help(name: &str) -> ExitCode {
    let mut cmd = Cli::command();
    // Build the tree so the subcommand knows its bin name for the usage line.
    cmd.build();
    if let Some(sub) = cmd.find_subcommand_mut(name) {
        let _ = sub.print_help();
    }
    ExitCode::from(2)
}

#[cfg(test)]
mod tests {
    use super::parse_num;

    #[test]
    fn parse_num_suffixes() {
        assert_eq!(parse_num("0"), Ok(0));
        assert_eq!(parse_num("250"), Ok(250));
        assert_eq!(parse_num("100k"), Ok(100_000));
        assert_eq!(parse_num("1M"), Ok(1_000_000));
        assert_eq!(parse_num("1.5k"), Ok(1500));
        assert_eq!(parse_num("2g"), Ok(2_000_000_000));
        assert!(parse_num("k").is_err() && parse_num("1x").is_err());
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    // One buffered stdout for every subcommand, flushed once below, so a failed
    // final write (full disk, closed pipe) is reported rather than lost when the
    // `BufWriter` drops.
    let stdout = std::io::stdout();
    let out = &mut std::io::BufWriter::new(stdout.lock());
    let result = match cli.command {
        Command::Extract {
            input,
            paf,
            reference,
            alt,
            source,
            min_clip,
            trim,
            min_gap,
            min_mapq,
            min_aln_len,
            flank,
            max_eseq,
            min_equal,
        } => {
            // With no input file, show the subcommand help rather than blocking
            // on stdin. Only an explicit "-" reads from stdin.
            let Some(input) = input else {
                return print_subcommand_help("extract");
            };
            let opts = ExtractOpts {
                input,
                paf,
                reference,
                alt,
                source,
                min_clip,
                trim,
                min_gap,
                min_mapq,
                min_aln_len,
                flank,
                max_eseq,
                min_equal,
            };
            extract::run(&opts, out)
        }
        Command::Getfa { min_equal, input } => {
            let Some(input) = input else {
                return print_subcommand_help("getfa");
            };
            getfa::run(&input, min_equal, out)
        }
        Command::Fltbnd {
            margin,
            source,
            min_equal,
            extract_out,
            rb3_paf,
        } => {
            let (Some(extract_out), Some(rb3_paf)) = (extract_out, rb3_paf) else {
                return print_subcommand_help("fltbnd");
            };
            let src = source.as_deref();
            fltbnd::run(&extract_out, &rb3_paf, margin, src, min_equal, out)
        }
        Command::Merge {
            inputs,
            min_cnt,
            min_cnt_strand,
            win_size,
            max_allele,
            max_check,
            min_equal,
            min_mapq,
            min_mapq_min,
            merge_input,
            min_cnt_in,
            min_cnt_strand_in,
        } => {
            if inputs.is_empty() {
                return print_subcommand_help("merge");
            }
            let opts = MergeOpts {
                inputs,
                min_cnt,
                min_cnt_strand,
                win_size,
                max_allele,
                max_check,
                min_equal,
                min_mapq,
                min_mapq_min,
                merge_input,
                min_cnt_in,
                min_cnt_strand_in,
            };
            merge::run(&opts, out)
        }
        Command::Fltreg { input, bed, margin } => {
            let (Some(input), Some(bed)) = (input, bed) else {
                return print_subcommand_help("fltreg");
            };
            fltreg::run(&input, &bed, margin, out)
        }
        Command::Fltcnt {
            input,
            min_src,
            min_src_strand,
            max_rest,
            min_dist,
            inv_dist,
            max_ins,
            src,
            rest,
        } => {
            let Some(input) = input else {
                return print_subcommand_help("fltcnt");
            };
            let opts = FltcntOpts {
                input,
                min_src,
                min_src_strand,
                max_rest,
                min_dist,
                inv_dist,
                max_ins,
                src,
                rest,
            };
            fltcnt::run(&opts, out)
        }
    };

    if let Err(e) = result.and_then(|()| out.flush()) {
        eprintln!("badclip: {e}");
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
