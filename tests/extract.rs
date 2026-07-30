//! End-to-end tests for `badclip extract` against the PAF/MSV fixtures.
//!
//! The provided `.msv` files carry a trailing INFO column produced by
//! minisv.js that differs from badclip's own INFO tags. So the expected output
//! is each `.msv` line truncated to its first 8 columns, with badclip's INFO
//! column (`aln_len=...;qlen=...;mapq=...`) appended per line.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn test_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test")
}

/// Expected output = each `.msv` line's first 8 columns, plus badclip's INFO
/// column `info[i]` (the entire 9th field) for output line `i`.
fn expected(name: &str, info: &[&str]) -> String {
    let msv = std::fs::read_to_string(test_dir().join(format!("{name}.msv"))).unwrap();
    let mut out = String::new();
    for (i, line) in msv.lines().enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        out.push_str(&cols[..8.min(cols.len())].join("\t"));
        out.push_str(&format!("\t{}\n", info[i]));
    }
    out
}

/// Run `badclip extract <paf>` and return its stdout.
fn run_extract_file(name: &str) -> String {
    run_extract_file_args(name, &[])
}

/// Run `badclip extract [args...] <paf>` and return its stdout.
fn run_extract_file_args(name: &str, args: &[&str]) -> String {
    let paf = test_dir().join(format!("{name}.paf"));
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .args(args)
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
    // Single hit (alen=14732); left clip (qs=141) then right clip (qlen-qe=292).
    assert_eq!(
        run_extract_file("clip01"),
        expected(
            "clip01",
            &[
                "aln_len=14732,0;qlen=14925,0,141;mapq=1,0",
                "aln_len=14732,0;qlen=14774,0,292;mapq=1,0",
            ]
        )
    );
}

#[test]
fn join01_matches() {
    // Flipped: chr22 (y1) first, chrY (y0) second.
    assert_eq!(
        run_extract_file("join01"),
        expected(
            "join01",
            &["aln_len=16881,2475;qlen=16848,15,2474;mapq=60,60"]
        )
    );
}

#[test]
fn join02_matches() {
    // Flipped: chr1 (y1) first, chr21 (y0) second.
    assert_eq!(
        run_extract_file("join02"),
        expected(
            "join02",
            &["aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60"]
        )
    );
}

#[test]
fn join03_matches() {
    // No flip: chr13 (y0) first, chr2 (y1) second.
    assert_eq!(
        run_extract_file("join03"),
        expected(
            "join03",
            &["aln_len=29505,5668;qlen=29436,1,5661;mapq=60,60"]
        )
    );
}

#[test]
fn min_aln_len_drops_all() {
    // -a above every hit's alignment length -> nothing survives.
    assert_eq!(run_extract_file_args("join01", &["-a", "100000"]), "");
}

#[test]
fn min_aln_len_drops_one_hit() {
    // -a 3000 drops join01's chrY hit (alen=2475); the surviving chr22 hit is
    // then a lone alignment, so only its left clip (qs=2489) is emitted.
    assert_eq!(
        run_extract_file_args("join01", &["-a", "3000"]),
        "chr22\t22131975\t<.\t.\t.\tm84039_230117_233243_s1/257233489/ccs\t60\t+\
         \taln_len=16881,0;qlen=16848,0,2489;mapq=60,0\n"
    );
}

#[test]
fn stdin_matches_file() {
    let paf = std::fs::read(test_dir().join("join02.paf")).unwrap();
    assert_eq!(
        run_extract_stdin(&paf),
        expected(
            "join02",
            &["aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60"]
        )
    );
}

#[test]
fn no_input_shows_help() {
    // `badclip extract` with no file must print help and exit 2, not block on
    // stdin. stdin is /dev/null so a regression would EOF rather than hang.
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .stdin(Stdio::null())
        .output()
        .expect("failed to run badclip");
    assert_eq!(output.status.code(), Some(2));
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(
        help.contains("Usage: badclip extract"),
        "expected help output, got: {help}"
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
    assert_eq!(
        run_extract_stdin(&gz),
        expected(
            "join01",
            &["aln_len=16881,2475;qlen=16848,15,2474;mapq=60,60"]
        )
    );
}
