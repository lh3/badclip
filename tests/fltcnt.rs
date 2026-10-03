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
    // -l -1 disables the distance/clip filter. Default -r 0: no --rest
    // (`s.raw`) reads allowed, so only the pure-`s.flt` tra passes -c 0.
    assert_eq!(kept(INPUT, &["-l", "-1", "-c", "0"]), ["tra"]);
    assert_eq!(kept(INPUT, &["-l", "-1"]), ["tra"]);
    // -r -1 lifts the cap. Default -c: >= 3 filtered reads.
    assert_eq!(kept(INPUT, &["-l", "-1", "-r", "-1"]), ["clip", "near", "far", "tra"]);
    assert_eq!(kept(INPUT, &["-l", "-1", "-r", "-1", "-c", "5"]), ["clip", "near", "far"]);
    assert_eq!(kept(INPUT, &["-l", "-1", "-r", "-1", "-c", "0"]).len(), 5);
    // -s: each filtered strand.
    assert_eq!(kept(INPUT, &["-l", "-1", "-r", "-1", "-s", "1"]), ["clip", "near"]);
    // -r: at most this many raw reads.
    assert_eq!(kept(INPUT, &["-l", "-1", "-c", "0", "-r", "1"]), ["far", "tra"]);
    assert_eq!(kept(INPUT, &["-l", "-1", "-c", "0", "-r", "3"]), ["clip", "near", "far", "tra"]);
}

#[test]
fn fltcnt_min_dist() {
    // Default -l 0 drops clips only.
    assert_eq!(kept(INPUT, &["-c", "0", "-r", "-1"]), ["near", "far", "tra", "raw"]);
    // Positive -l also drops close same-contig calls; inter-contig joins pass.
    assert_eq!(kept(INPUT, &["-c", "0", "-r", "-1", "-l", "100"]), ["far", "tra", "raw"]);
}

#[test]
fn fltcnt_explicit_sources() {
    let input = "\
chr1\t100\t>>\tchr1\t5000\tx\t6\t+\tavg_mapq=60,60;count=A:1,1|B:3,3|C:2,0
";
    // No *.flt source: everything is --rest, so -c fails, with a warning.
    let out = run_fltcnt(input, &[]);
    assert!(out.status.success() && out.stdout.is_empty());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--src"));
    // --src B (6 reads); --rest defaults to A+C (4 reads), over the -r 0 cap.
    assert!(kept(input, &["--src", "B"]).is_empty());
    assert!(kept(input, &["--src", "B", "-r", "3"]).is_empty());
    assert_eq!(kept(input, &["--src", "B", "-r", "4"]), ["x"]);
    assert_eq!(kept(input, &["--src", "B", "-r", "-1"]), ["x"]);
    // Comma lists: --src A,C has 4 reads (fails -c 5), and C has 0 on reverse.
    assert_eq!(kept(input, &["--src", "A,C", "-r", "-1"]), ["x"]);
    assert!(kept(input, &["--src", "A,C", "-c", "5", "-r", "-1"]).is_empty());
    assert!(kept(input, &["--src", "A,C", "-c", "4", "-s", "2", "-r", "-1"]).is_empty());
    // All sources in --src: nothing left for --rest.
    assert_eq!(kept(input, &["--src", "A,B,C"]), ["x"]);
    // Explicit --rest: only A counts (2 reads); C is in neither group.
    assert!(kept(input, &["--src", "B", "--rest", "A", "-r", "1"]).is_empty());
    assert_eq!(kept(input, &["--src", "B", "--rest", "A", "-r", "2"]), ["x"]);
}

#[test]
fn fltcnt_no_source_aborts() {
    // `merge -m` output has no source labels: neither group matches.
    let input = "chr1\t100\t>>\tchr1\t5000\tx\t6\t+\tavg_mapq=60,60;count=3,3\n";
    let out = run_fltcnt(input, &[]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--src"));
    // With explicit lists, a line with only other sources also aborts.
    let input = "chr1\t100\t>>\tchr1\t5000\tx\t6\t+\tavg_mapq=60,60;count=Z:3,3\n";
    assert!(!run_fltcnt(input, &["--src", "B", "--rest", "A"]).status.success());
}
