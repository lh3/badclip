//! badclip — extract breakends and structural-variant signals from long-read
//! alignments.

mod aln;
mod extract;
mod fltreg;
mod flteseq;
mod geteseq;
#[allow(dead_code)]
mod iitree;
mod io;
mod merge;
mod paf;

use std::process::ExitCode;

use clap::{CommandFactory, Parser, Subcommand};

use extract::ExtractOpts;
use merge::MergeOpts;

#[derive(Parser)]
#[command(name = "badclip", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Extract breakends from read alignments in BAM (or PAF) format.
    Extract {
        /// Input alignments: BAM by default, or PAF with --paf. Use "-" for stdin.
        input: Option<String>,

        /// Read PAF (optionally gzip'd) instead of BAM.
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

        /// Minimum clip length to report a clip breakend.
        #[arg(short = 'c', long = "min-clip", default_value_t = 50)]
        min_clip: i64,

        /// Drop hits whose alignment block length (PAF col 11) is below this.
        #[arg(short = 'a', long = "min-aln-len", default_value_t = 0)]
        min_aln_len: i64,

        /// Flanking read sequence extracted on each side of a breakend (BAM only).
        #[arg(short = 'f', long = "flank", default_value_t = 250)]
        flank: i64,

        /// Maximum extracted sequence length; longer windows omit eseq (BAM only).
        #[arg(short = 'e', long = "max-eseq", default_value_t = 1000)]
        max_eseq: i64,

        /// Drop breakends whose eseq quality (`equal`) is below this (0 = keep all).
        #[arg(short = 'Q', long = "min-equal", default_value_t = 0)]
        min_equal: i64,

        /// Drop output lines whose col-7 mapq is below this (post-filter).
        #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
        min_mapq: i64,
    },

    /// Convert `extract` output into a FASTA of extracted breakend sequences.
    Geteseq {
        /// Drop records whose eseq quality (`equal` tag) is below this.
        #[arg(short = 'Q', long = "min-equal", default_value_t = 20)]
        min_equal: i64,

        /// Input `extract` output. Omit or use "-" to read from stdin.
        input: Option<String>,
    },

    /// Filter `extract` breakends whose eseq is spanned by a pangenome alignment.
    Flteseq {
        /// Margin extending the protected junction interval on each side.
        #[arg(short = 'l', long = "margin", default_value_t = 50)]
        margin: i64,

        /// Print all input lines, rewriting survivors' `source=` tag to this.
        #[arg(short = 's', long = "source")]
        source: Option<String>,

        /// Drop lines whose eseq quality (`equal` tag) is below this.
        #[arg(short = 'Q', long = "min-equal", default_value_t = 20)]
        min_equal: i64,

        /// `extract` output (gzip ok; "-" for stdin).
        extract_out: Option<String>,

        /// ropebwt3 `sw` PAF on the `geteseq` FASTA (gzip ok).
        rb3_paf: Option<String>,
    },

    /// Merge per-read `extract` breakends into consensus SV calls.
    Merge {
        /// `extract` output (gzip ok; "-" for stdin).
        input: Option<String>,

        /// Minimum read count to emit a call.
        #[arg(short = 'c', long = "min-cnt", default_value_t = 3)]
        min_cnt: i64,

        /// Minimum read count on each strand.
        #[arg(short = 's', long = "min-cnt-strand", default_value_t = 1)]
        min_cnt_strand: i64,

        /// Clustering window size (bp).
        #[arg(short = 'w', long = "win-size", default_value_t = 100)]
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

        /// Drop input breakends whose col-7 mapq is below this.
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
    },
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

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Extract {
            input,
            paf,
            reference,
            alt,
            source,
            min_clip,
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
            extract::run(&ExtractOpts {
                input,
                paf,
                reference,
                alt,
                source,
                min_clip,
                min_mapq,
                min_aln_len,
                flank,
                max_eseq,
                min_equal,
            })
        }
        Command::Geteseq { min_equal, input } => {
            let Some(input) = input else {
                return print_subcommand_help("geteseq");
            };
            geteseq::run(&input, min_equal)
        }
        Command::Flteseq {
            margin,
            source,
            min_equal,
            extract_out,
            rb3_paf,
        } => {
            let (Some(extract_out), Some(rb3_paf)) = (extract_out, rb3_paf) else {
                return print_subcommand_help("flteseq");
            };
            flteseq::run(&extract_out, &rb3_paf, margin, source.as_deref(), min_equal)
        }
        Command::Merge {
            input,
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
            let Some(input) = input else {
                return print_subcommand_help("merge");
            };
            merge::run(&MergeOpts {
                input,
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
            })
        }
        Command::Fltreg { input, bed } => {
            let (Some(input), Some(bed)) = (input, bed) else {
                return print_subcommand_help("fltreg");
            };
            fltreg::run(&input, &bed)
        }
    };

    if let Err(e) = result {
        eprintln!("badclip: {e}");
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
