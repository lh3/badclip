//! End-to-end tests for `badclip fltreg` against `test/fltreg01.bed`.
//!
//! The fixture has overlapping regions `chr1 100-200` and `chr1 150-250` (which
//! must merge to `chr1 [100,250)`), `chr1 [5000,6000)`, and `chr2 [50,60)`.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn test_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test")
}

/// Run `badclip fltreg - <bed>` feeding `input` on stdin; return stdout.
fn run_fltreg(input: &str) -> String {
    let bed = test_dir().join("fltreg01.bed");
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("fltreg")
        .arg("-")
        .arg(&bed)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip fltreg");
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
fn fltreg_filters_either_breakend() {
    // A line is dropped if either endpoint falls in a region.
    let input = "\
chr1\t150\t>>\tchr1\t9000\tep1_in\t60\t+\tsrc
chr1\t9000\t>>\tchr1\t5500\tep2_in\t60\t+\tsrc
chr1\t240\t>.\t.\t.\tclip_in_merged\t60\t+\tsrc
chr2\t55\t>>\tchr2\t9000\tctg2_in\t60\t+\tsrc
chr1\t9000\t>>\tchr1\t9500\tboth_out\t60\t+\tsrc
chr3\t150\t>>\tchr3\t150\tother_ctg\t60\t+\tsrc
";
    let got = run_fltreg(input);
    // ep1_in (150 in [100,250)), ep2_in (5500 in [5000,6000)), clip_in_merged
    // (240 in the merged [100,250)), and ctg2_in (55 in [50,60)) are dropped.
    for dropped in ["ep1_in", "ep2_in", "clip_in_merged", "ctg2_in"] {
        assert!(!got.contains(dropped), "{dropped} should be dropped:\n{got}");
    }
    // both_out and other_ctg (contig not in BED) survive.
    assert!(got.contains("both_out") && got.contains("other_ctg"), "survivors wrong:\n{got}");
    assert_eq!(got.lines().count(), 2, "expected two survivors:\n{got}");
}

#[test]
fn fltreg_half_open_boundary() {
    // [start,end): p==start is inside (dropped), p==end is outside (kept).
    let input = "\
chr1\t100\t>.\t.\t.\tat_start\t60\t+\tsrc
chr1\t250\t>.\t.\t.\tat_end\t60\t+\tsrc
";
    let got = run_fltreg(input);
    assert!(!got.contains("at_start"), "p==start should be dropped:\n{got}");
    assert!(got.contains("at_end"), "p==end (merged) should be kept:\n{got}");
    assert_eq!(got.lines().count(), 1, "only at_end should survive:\n{got}");
}

#[test]
fn fltreg_survivors_are_verbatim() {
    // A surviving line is printed byte-for-byte unchanged.
    let line = "chr1\t9000\t>>\tchr1\t9500\tr\t60\t+\tsource=foo;idx=0;mapq=60,60";
    let got = run_fltreg(&format!("{line}\n"));
    assert_eq!(got, format!("{line}\n"));
}

#[test]
fn fltreg_no_input_shows_help() {
    // No file args -> print help and exit 2, not block on stdin.
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("fltreg")
        .stdin(Stdio::null())
        .output()
        .expect("failed to run badclip fltreg");
    assert_eq!(output.status.code(), Some(2));
    assert!(
        String::from_utf8(output.stdout).unwrap().contains("Usage: badclip fltreg"),
        "expected usage message"
    );
}
