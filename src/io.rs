//! Input helpers: seamlessly read from stdin or a file, transparently
//! decompressing gzip when the input begins with the gzip magic bytes.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::sync::atomic::{AtomicBool, Ordering};

use flate2::read::MultiGzDecoder;

/// Open `path` as a buffered reader.
///
/// - `"-"` reads from standard input (at most one input per run; [`open_raw`]).
/// - Any other value is treated as a file path.
///
/// In both cases the first two bytes are peeked (without consuming them); if
/// they are the gzip magic `1f 8b` the stream is wrapped in a
/// [`MultiGzDecoder`], so gzip'd and plain inputs are handled identically.
pub fn open_reader(path: &str) -> io::Result<Box<dyn BufRead>> {
    let mut reader = BufReader::new(open_raw(path)?);
    // `fill_buf` peeks the buffered bytes without consuming them, so the
    // magic-byte check does not disturb the stream that is later decoded.
    let is_gzip = {
        let head = reader.fill_buf()?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };

    if is_gzip {
        // MultiGzDecoder handles concatenated members (e.g. bgzf blocks).
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

/// Whether some input has already claimed stdin (`"-"`).
static STDIN_TAKEN: AtomicBool = AtomicBool::new(false);

/// Open `path` unbuffered: `"-"` is stdin (locked), anything else a file. Only
/// one input per run can be `"-"`; a second is an error, since the first would
/// have consumed stdin (or still hold its lock, and hang).
pub fn open_raw(path: &str) -> io::Result<Box<dyn Read>> {
    if path != "-" {
        return Ok(Box::new(File::open(path)?));
    }
    if STDIN_TAKEN.swap(true, Ordering::Relaxed) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only one input can be `-` (stdin)",
        ));
    }
    Ok(Box::new(io::stdin().lock()))
}

/// `path` for messages: `"-"` reads as `stdin`.
pub fn display_path(path: &str) -> &str {
    if path == "-" { "stdin" } else { path }
}

/// An `InvalidData` error at line `lineno` (1-based) of input `path`.
pub fn line_error(path: &str, lineno: usize, msg: impl std::fmt::Display) -> io::Error {
    let path = display_path(path);
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("{path}:{lineno}: {msg}"),
    )
}
