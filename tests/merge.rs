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
    // Default thresholds (-c 3 -s 1): the >> join cluster (A), the clip cluster
    // (C), and the inversion-pair join cluster (D, with count_fr/count_rf) all
    // pass; the all-`+` B cluster fails the per-strand filter (0 reads on `-`),
    // and the singleton readSolo and the 2-read E cluster fail min_cnt.
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
fn merge_min_mapq_filters() {
    // Three clips clustering at chr1:1000; qlo has col-7 mapq 5. -q drops input
    // breakends whose col-7 mapq is below the threshold before clustering.
    let input = "\
chr1\t1000\t>.\t.\t.\tqhi\t60\t+\tsource=foo
chr1\t1000\t>.\t.\t.\tqhi2\t60\t-\tsource=foo
chr1\t1000\t>.\t.\t.\tqlo\t5\t+\tsource=foo
";
    // -q 20 drops qlo (mapq 5); qhi/qhi2 remain.
    let got = run_merge_stdin(input, &["-c", "1", "-s", "1", "-q", "20"]);
    assert!(
        got.contains("qhi") && got.contains("qhi2") && !got.contains("qlo"),
        "qlo (mapq 5) should be filtered by -q 20:\n{got}"
    );
    // Default -q 0 keeps everything, so qlo appears.
    let all = run_merge_stdin(input, &["-c", "1", "-s", "1"]);
    assert!(all.contains("qlo"), "qlo should survive default -q 0:\n{all}");
}

#[test]
fn merge_min_mapq_min_filters() {
    // Two joins and one clip clustering at chr1:1000->chr2:2000. jlo's `mapq=`
    // pair has a smaller value of 3; jhi's is 40. -p drops join breakends whose
    // *smaller* per-hit mapq is below the threshold, but never a clip (mapq2=0).
    let input = "\
chr1\t1000\t>>\tchr2\t2000\tjhi\t60\t+\tsource=foo;mapq=60,40
chr1\t1000\t>>\tchr2\t2000\tjhi2\t60\t-\tsource=foo;mapq=40,60
chr1\t1000\t>>\tchr2\t2000\tjlo\t60\t+\tsource=foo;mapq=60,3
chr1\t3000\t>.\t.\t.\tcliplo\t60\t+\tsource=foo;mapq=3,0
chr1\t3000\t>.\t.\t.\tcliplo2\t60\t-\tsource=foo;mapq=3,0
";
    // -p 20 drops jlo (min pair 3); jhi/jhi2 remain and the clip is exempt.
    let got = run_merge_stdin(input, &["-c", "1", "-s", "1", "-p", "20"]);
    assert!(
        got.contains("jhi") && got.contains("jhi2") && !got.contains("jlo"),
        "jlo (min pair mapq 3) should be filtered by -p 20:\n{got}"
    );
    assert!(
        got.contains("cliplo"),
        "a clip (mapq2=0) must be exempt from -p:\n{got}"
    );
    // -p 0 keeps every join, so jlo appears.
    let all = run_merge_stdin(input, &["-c", "1", "-s", "1", "-p", "0"]);
    assert!(all.contains("jlo"), "jlo should survive -p 0:\n{all}");
}

#[test]
fn merge_avg_mapq_per_end() {
    // avg_mapq=q1,q2 averages the `mapq=` pair per output endpoint: q1 on the
    // ctg1:pos1 side, q2 on the ctg2:pos2 side. The cluster mixes ">< " and "<>"
    // (an inversion pair) with asymmetric per-hit mapqs — since the INFO pair is
    // already endpoint-aligned, the "<>" member must NOT swap sides.
    let input = "\
chr1\t1000\t><\tchr1\t5000\tri1\t60\t+\tsource=foo;mapq=60,20
chr1\t1005\t<>\tchr1\t5005\tri2\t62\t-\tsource=foo;mapq=62,22
chr1\t1010\t><\tchr1\t4995\tri3\t58\t+\tsource=foo;mapq=58,18
";
    // q1 = mean(60,62,58) = 60; q2 = mean(20,22,18) = 20.
    let got = run_merge_stdin(input, &["-p", "0"]);
    assert!(
        got.contains("avg_mapq=60,20"),
        "per-end avg_mapq should be 60,20 (endpoint-aligned across the inversion pair):\n{got}"
    );

    // A clip cluster reports q2=0 (no second end).
    let clips = "\
chr1\t2000\t>.\t.\t.\trc1\t40\t+\tsource=foo;mapq=40,0
chr1\t2005\t>.\t.\t.\trc2\t20\t-\tsource=foo;mapq=20,0
chr1\t2010\t>.\t.\t.\trc3\t60\t+\tsource=foo;mapq=60,0
";
    let cg = run_merge_stdin(clips, &[]);
    assert!(
        cg.contains("avg_mapq=40,0"),
        "clip cluster avg_mapq should be 40,0:\n{cg}"
    );
}

#[test]
fn merge_sample_merge() {
    // -m: input is combined `merge` output (one line per sample). Two samples'
    // join calls cluster at chr1; two inversion-pair calls cluster at chr2.
    // Counts come from the `count=` tag (summed over sources), mapq from
    // `avg_mapq=`; output has count=F,R, no per-source breakdown, no reads=.
    let input = "\
chr1\t1000\t>>\tchr1\t5000\t.\t6\t+\tavg_mapq=60,50;count=sA:2,4
chr1\t1010\t>>\tchr1\t5008\t.\t5\t-\tavg_mapq=40,30;count=sB:3,2
chr2\t2000\t><\tchr2\t2800\t.\t6\t+\tavg_mapq=60,60;count=sA:3,3;count_fr=6;count_rf=0
chr2\t2005\t<>\tchr2\t2795\t.\t6\t-\tavg_mapq=60,60;count=sB:2,4;count_fr=2;count_rf=4
";
    let got = run_merge_stdin(input, &["-m"]);
    // chr1: rep = 2nd member (pos 1010); q1=mean(60,40)=50, q2=mean(50,30)=40;
    // fwd=2+3=5, rev=4+2=6; col-7 = 2 distinct sources (sA, sB).
    assert!(
        got.contains("chr1\t1010\t>>\tchr1\t5008\t.\t2\t-\tavg_mapq=50,40;count=5,6"),
        "chr1 sample-merge line wrong:\n{got}"
    );
    // chr2 inversion: count_fr/count_rf summed from input (6+2, 0+4); fr*rf != 0
    // so no foldback despite ctg1==ctg2; col-7 = 2 sources.
    assert!(
        got.contains("chr2\t2005\t<>\tchr2\t2795\t.\t2\t-\tavg_mapq=60,60;count=5,7;count_fr=8;count_rf=4"),
        "chr2 inversion sample-merge line wrong:\n{got}"
    );
    assert!(!got.contains("foldback"), "unexpected foldback (fr*rf != 0):\n{got}");
    assert!(!got.contains("reads="), "-m output must not carry a reads= tag:\n{got}");
    assert_eq!(got.lines().count(), 2, "expected two merged clusters:\n{got}");
}

#[test]
fn merge_sample_merge_per_line_filters() {
    // -C/-S drop under-supported *input* lines (before clustering), -m only.
    // Three clips at chr3: sA total 1 (< -C 3), sB rev 0 (< -S 1), sC ok.
    let input = "\
chr3\t3000\t>.\t.\t.\t.\t5\t+\tavg_mapq=40,0;count=sA:1,0
chr3\t3000\t>.\t.\t.\t.\t5\t-\tavg_mapq=40,0;count=sB:5,0
chr3\t3005\t>.\t.\t.\t.\t5\t+\tavg_mapq=40,0;count=sC:2,2
";
    // Defaults -C 3 -S 1: only sC survives parsing -> one cluster, count=2,2.
    let got = run_merge_stdin(input, &["-m"]);
    assert_eq!(got.lines().count(), 1, "only sC should survive -C/-S:\n{got}");
    assert!(got.contains("count=2,2"), "surviving cluster wrong:\n{got}");
    // -C 0 -S 0 keeps all three; they cluster into one: fwd=1+5+2=8, rev=0+0+2=2,
    // col-7 = 3 distinct sources (sA, sB, sC).
    let all = run_merge_stdin(input, &["-m", "-C", "0", "-S", "0"]);
    assert!(
        all.contains("\t3\t-\tavg_mapq=40,0;count=8,2"),
        "relaxed -C/-S should keep all three (col-7 = 3 sources):\n{all}"
    );
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
