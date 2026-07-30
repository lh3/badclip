//! badclip — extract breakends and structural-variant signals from long-read
//! alignments.

mod bam;
mod extract;
mod geteseq;
mod io;
mod paf;

use std::process::ExitCode;

use clap::{CommandFactory, Parser, Subcommand};

use extract::ExtractOpts;

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

        /// Minimum clip length to report a clip breakend.
        #[arg(short = 'c', long = "min-clip", default_value_t = 100)]
        min_clip: i64,

        /// Drop hits with mapping quality below this value.
        #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
        min_mapq: i64,

        /// Drop hits whose alignment block length (PAF col 11) is below this.
        #[arg(short = 'a', long = "min-aln-len", default_value_t = 0)]
        min_aln_len: i64,

        /// Flanking read sequence extracted on each side of a breakend (BAM only).
        #[arg(short = 'f', long = "flank", default_value_t = 250)]
        flank: i64,

        /// Maximum extracted sequence length; longer windows omit eseq (BAM only).
        #[arg(short = 'e', long = "max-eseq", default_value_t = 1000)]
        max_eseq: i64,
    },

    /// Convert `extract` output into a FASTA of extracted breakend sequences.
    Geteseq {
        /// Input `extract` output. Omit or use "-" to read from stdin.
        input: Option<String>,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Extract {
            input,
            paf,
            min_clip,
            min_mapq,
            min_aln_len,
            flank,
            max_eseq,
        } => {
            // With no input file, show the subcommand help rather than blocking
            // on stdin. Only an explicit "-" reads from stdin.
            let Some(input) = input else {
                let mut cmd = Cli::command();
                // Build the tree so subcommands know their bin name ("badclip
                // extract") for the rendered usage line.
                cmd.build();
                if let Some(sub) = cmd.find_subcommand_mut("extract") {
                    let _ = sub.print_help();
                }
                return ExitCode::from(2);
            };
            extract::run(&ExtractOpts {
                input,
                paf,
                min_clip,
                min_mapq,
                min_aln_len,
                flank,
                max_eseq,
            })
        }
        Command::Geteseq { input } => geteseq::run(input.as_deref()),
    };

    if let Err(e) = result {
        eprintln!("badclip: {e}");
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
