//! End-to-end tests for `badclip fltcnt` on inline `merge`-shaped input.

use std::io::Write;
use std::process::{Command, Output, Stdio};

/// Run `badclip fltcnt [args..] -` feeding `input` on stdin.
fn run_fltcnt(input: &str, args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("fltcnt")
        .args(args)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip fltcnt");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Names (col 6 slot reused as a label) of the lines that survive.
fn kept(input: &str, args: &[&str]) -> Vec<String> {
    let out = run_fltcnt(input, args);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').nth(5).unwrap().to_string())
        .collect()
}

const INPUT: &str = "\
chr1\t100\t>.\t.\t.\tclip\t9\t+\tavg_mapq=60,0;count=s.flt:3,3|s.raw:2,1
chr1\t100\t>>\tchr1\t150\tnear\t9\t+\tavg_mapq=60,60;count=s.flt:3,3|s.raw:2,1
chr1\t100\t>>\tchr1\t5000\tfar\t6\t+\tavg_mapq=60,60;count=s.flt:5,0|s.raw:1,0
chr1\t100\t>>\tchr2\t50\ttra\t4\t+\tavg_mapq=60,60;count=s.flt:4,0
chr1\t100\t>>\tchr1\t9000\traw\t8\t+\tavg_mapq=60,60;count=s.raw:4,4
";

#[test]
fn fltcnt_counts() {
    // -l -1 disables the distance/clip filter. Default -c: >= 5 filtered reads.
    assert_eq!(kept(INPUT, &["-l", "-1"]), ["clip", "near", "far"]);
    assert_eq!(kept(INPUT, &["-l", "-1", "-c", "0"]).len(), 5);
    // -s: each filtered strand.
    assert_eq!(kept(INPUT, &["-l", "-1", "-s", "1"]), ["clip", "near"]);
    // -r: raw reads.
    assert_eq!(kept(INPUT, &["-l", "-1", "-c", "0", "-r", "3"]), ["clip", "near", "raw"]);
}

#[test]
fn fltcnt_min_dist() {
    // Default -l 0 drops clips only.
    assert_eq!(kept(INPUT, &["-c", "0"]), ["near", "far", "tra", "raw"]);
    // Positive -l also drops close same-contig calls; inter-contig joins pass.
    assert_eq!(kept(INPUT, &["-c", "0", "-l", "100"]), ["far", "tra", "raw"]);
}

#[test]
fn fltcnt_explicit_sources() {
    let input = "chr1\t100\t>>\tchr1\t5000\tx\t6\t+\tavg_mapq=60,60;count=A:1,1|B:3,3\n";
    // Neither source matches the *.raw/*.flt convention: abort.
    let out = run_fltcnt(input, &[]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--src-raw"));
    // Explicit names: B is the filtered source (6 reads), A the raw one (2).
    assert_eq!(kept(input, &["--src-raw", "A", "--src-flt", "B"]), ["x"]);
    assert!(kept(input, &["--src-raw", "A", "--src-flt", "B", "-r", "3"]).is_empty());
}
