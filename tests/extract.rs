//! End-to-end tests for `badclip extract` against the PAF/MSV fixtures.
//!
//! The provided `.msv` files carry a trailing INFO column produced by
//! minisv.js that differs from badclip's own `aln_len` tag (both in value —
//! query-span vs. PAF col-11 block length — and ordering). So the expected
//! output is each `.msv` line truncated to its first 8 columns, with badclip's
//! `aln_len=...` INFO appended per line.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn test_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test")
}

/// Expected output = each `.msv` line's first 8 columns, plus `aln_len=<suffix>`
/// where `alns[i]` is the expected `aln_len` value for output line `i`.
fn expected(name: &str, alns: &[&str]) -> String {
    let msv = std::fs::read_to_string(test_dir().join(format!("{name}.msv"))).unwrap();
    let mut out = String::new();
    for (i, line) in msv.lines().enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        out.push_str(&cols[..8.min(cols.len())].join("\t"));
        out.push_str(&format!("\taln_len={}\n", alns[i]));
    }
    out
}

/// Run `badclip extract <paf>` and return its stdout.
fn run_extract_file(name: &str) -> String {
    let paf = test_dir().join(format!("{name}.paf"));
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg(&paf)
        .output()
        .expect("failed to run badclip");
    assert!(output.status.success(), "badclip exited with failure");
    String::from_utf8(output.stdout).unwrap()
}

/// Run `badclip extract -` feeding `bytes` on stdin.
fn run_extract_stdin(bytes: &[u8]) -> String {
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip");
    child.stdin.take().unwrap().write_all(bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "badclip exited with failure");
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn clip01_matches() {
    // Single hit (alen=14732); both clip lines carry len2=0.
    assert_eq!(
        run_extract_file("clip01"),
        expected("clip01", &["14732,0", "14732,0"])
    );
}

#[test]
fn join01_matches() {
    // Flipped: chr22 (y1, alen=16881) first, chrY (y0, alen=2475) second.
    assert_eq!(
        run_extract_file("join01"),
        expected("join01", &["16881,2475"])
    );
}

#[test]
fn join02_matches() {
    // Flipped: chr1 (y1, alen=16879) first, chr21 (y0, alen=26284) second.
    assert_eq!(
        run_extract_file("join02"),
        expected("join02", &["16879,26284"])
    );
}

#[test]
fn join03_matches() {
    // No flip: chr13 (y0, alen=29505) first, chr2 (y1, alen=5668) second.
    assert_eq!(
        run_extract_file("join03"),
        expected("join03", &["29505,5668"])
    );
}

#[test]
fn stdin_matches_file() {
    let paf = std::fs::read(test_dir().join("join02.paf")).unwrap();
    assert_eq!(
        run_extract_stdin(&paf),
        expected("join02", &["16879,26284"])
    );
}

#[test]
fn gzip_is_detected() {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    let paf = std::fs::read(test_dir().join("join01.paf")).unwrap();
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(&paf).unwrap();
    let gz = enc.finish().unwrap();
    assert_eq!(run_extract_stdin(&gz), expected("join01", &["16881,2475"]));
}
