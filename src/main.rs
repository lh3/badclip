//! badclip — extract breakends and structural-variant signals from long-read
//! alignments.

mod extract;
mod io;
mod paf;

use std::process::ExitCode;

use clap::{Parser, Subcommand};

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
        /// Input PAF (optionally gzip'd). Omit or use "-" to read from stdin.
        input: Option<String>,

        /// Minimum clip length to report a clip breakend.
        #[arg(short = 'c', long = "min-clip", default_value_t = 100)]
        min_clip: i64,

        /// Drop hits with mapping quality below this value.
        #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
        min_mapq: i64,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Extract {
            input,
            min_clip,
            min_mapq,
        } => extract::run(&ExtractOpts {
            input,
            min_clip,
            min_mapq,
        }),
    };

    if let Err(e) = result {
        eprintln!("badclip: {e}");
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
