//! badclip — extract breakends and structural-variant signals from long-read
//! alignments.

mod extract;
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
    /// Extract breakends from read alignments in PAF format.
    Extract {
        /// Input PAF (optionally gzip'd). Use "-" to read from stdin.
        input: Option<String>,

        /// Minimum clip length to report a clip breakend.
        #[arg(short = 'c', long = "min-clip", default_value_t = 100)]
        min_clip: i64,

        /// Drop hits with mapping quality below this value.
        #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
        min_mapq: i64,

        /// Drop hits whose alignment block length (PAF col 11) is below this.
        #[arg(short = 'a', long = "min-aln-len", default_value_t = 0)]
        min_aln_len: i64,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Extract {
            input,
            min_clip,
            min_mapq,
            min_aln_len,
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
                min_clip,
                min_mapq,
                min_aln_len,
            })
        }
    };

    if let Err(e) = result {
        eprintln!("badclip: {e}");
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
