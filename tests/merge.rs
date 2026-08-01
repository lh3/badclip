//! End-to-end tests for `badclip merge` against the hand-authored fixture.
//!
//! `test/merge01.clip` is a small, deliberately-unsorted `extract`-format input
//! (exercising the in-memory sort). `test/merge01.expected` is the golden output
//! under default thresholds, hand-derived from the documented format (NOT diffed
//! against `minisv.js merge`, which we intentionally diverge from).

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn test_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test")
}

/// Run `badclip merge [args...] <file>` and return its stdout.
fn run_merge(args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("merge")
        .args(args)
        .output()
        .expect("failed to run badclip merge");
    assert!(
        output.status.success(),
        "badclip merge failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn merge_basic() {
    // Default thresholds (-c 4 -s 2): the >> join cluster (A), the clip cluster
    // (C), and the inversion-pair join cluster (D, with count_fr/count_rf) all
    // pass; the all-`+` B cluster fails the per-strand filter, and the singletons
    // (readSolo, the 2-read E cluster) fail min_cnt.
    let got = run_merge(&[test_dir().join("merge01.clip").to_str().unwrap()]);
    let want = std::fs::read_to_string(test_dir().join("merge01.expected")).unwrap();
    assert_eq!(got, want);
}

#[test]
fn merge_thresholds_relax() {
    // Lowering both thresholds (-c 2 -s 1) lets the balanced 2-read E cluster
    // through in addition to the three default calls.
    let got = run_merge(&[
        "-c",
        "2",
        "-s",
        "1",
        test_dir().join("merge01.clip").to_str().unwrap(),
    ]);
    assert_eq!(got.lines().count(), 4);
    assert!(
        got.contains("reads=bar:readE2|foo:readE1"),
        "expected the 2-read E cluster (source-stratified reads) to survive -c 2 -s 1:\n{got}"
    );
}

#[test]
fn merge_reads_stdin() {
    // `-` reads from stdin, and gives the same result as the file path.
    let bytes = std::fs::read(test_dir().join("merge01.clip")).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("merge")
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip merge");
    child.stdin.take().unwrap().write_all(&bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let want = std::fs::read_to_string(test_dir().join("merge01.expected")).unwrap();
    assert_eq!(String::from_utf8(output.stdout).unwrap(), want);
}

/// Run `badclip merge [args...] -` feeding `input` on stdin.
fn run_merge_stdin(input: &str, args: &[&str]) -> String {
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("merge")
        .args(args)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip merge");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(out.status.success());
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn merge_min_equal_filters() {
    // Three clips clustering at chr1:1000; qlow has equal=5. -Q drops input
    // breakends whose `equal` is below the threshold before clustering.
    let input = "\
chr1\t1000\t>.\t.\t.\tqfwd\t60\t+\tsource=foo;equal=30
chr1\t1000\t>.\t.\t.\tqrev\t60\t-\tsource=foo;equal=30
chr1\t1000\t>.\t.\t.\tqlow\t60\t+\tsource=foo;equal=5
";
    // Default -Q 20: qlow (equal=5) is dropped before clustering.
    let got = run_merge_stdin(input, &["-c", "1", "-s", "1"]);
    assert!(
        got.contains("qfwd") && got.contains("qrev") && !got.contains("qlow"),
        "qlow (equal=5) should be filtered by default -Q 20:\n{got}"
    );
    // -Q 0 keeps everything, so qlow appears.
    let all = run_merge_stdin(input, &["-c", "1", "-s", "1", "-Q", "0"]);
    assert!(all.contains("qlow"), "qlow should survive -Q 0:\n{all}");
}

#[test]
fn merge_no_input_shows_help() {
    // No file -> print help and exit 2, not block on stdin (stdin is /dev/null).
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("merge")
        .stdin(Stdio::null())
        .output()
        .expect("failed to run badclip merge");
    assert_eq!(output.status.code(), Some(2));
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(
        help.contains("Usage: badclip merge"),
        "expected usage, got: {help}"
    );
}
