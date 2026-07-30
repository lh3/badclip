//! Input helpers: seamlessly read from stdin or a file, transparently
//! decompressing gzip when the input begins with the gzip magic bytes.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

use flate2::read::MultiGzDecoder;

/// Open `path` as a buffered reader.
///
/// - `"-"` reads from standard input.
/// - Any other value is treated as a file path.
///
/// In both cases the first two bytes are peeked (without consuming them); if
/// they are the gzip magic `1f 8b` the stream is wrapped in a
/// [`MultiGzDecoder`], so gzip'd and plain inputs are handled identically.
pub fn open_reader(path: &str) -> io::Result<Box<dyn BufRead>> {
    let raw: Box<dyn Read> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path)?)
    };

    let mut reader = BufReader::new(raw);
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
