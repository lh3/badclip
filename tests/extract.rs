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

/// Run `badclip extract --paf <paf>` and return its stdout.
fn run_extract_file(name: &str) -> String {
    run_extract_file_args(name, &[])
}

/// Run `badclip extract --paf [args...] <paf>` and return its stdout.
fn run_extract_file_args(name: &str, args: &[&str]) -> String {
    let paf = test_dir().join(format!("{name}.paf"));
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg("--paf")
        .args(args)
        .arg(&paf)
        .output()
        .expect("failed to run badclip");
    assert!(output.status.success(), "badclip exited with failure");
    String::from_utf8(output.stdout).unwrap()
}

/// Run `badclip extract <bam>` and return its stdout (BAM is the default input).
fn run_extract_bam(name: &str) -> String {
    run_extract_bam_args(name, &[])
}

/// Run `badclip extract [args...] <bam>` and return its stdout.
fn run_extract_bam_args(name: &str, args: &[&str]) -> String {
    let bam = test_dir().join(format!("{name}.bam"));
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .args(args)
        .arg(&bam)
        .output()
        .expect("failed to run badclip");
    assert!(
        output.status.success(),
        "badclip failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

/// Run `badclip extract --paf -` feeding `bytes` on stdin.
fn run_extract_stdin(bytes: &[u8]) -> String {
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg("--paf")
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
                "source=foo;idx=0;aln_len=14732,0;qlen=14925,0,141;mapq=1,0",
                "source=foo;idx=1;aln_len=14732,0;qlen=14774,0,292;mapq=1,0",
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
            &["source=foo;idx=0;aln_len=16881,2475;qlen=16848,15,2474;mapq=60,60"]
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
            &["source=foo;idx=0;aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60"]
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
            &["source=foo;idx=0;aln_len=29505,5668;qlen=29436,1,5661;mapq=60,60"]
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
         \tsource=foo;idx=0;aln_len=16881,0;qlen=16848,0,2489;mapq=60,0\n"
    );
}

#[test]
fn stdin_matches_file() {
    let paf = std::fs::read(test_dir().join("join02.paf")).unwrap();
    assert_eq!(
        run_extract_stdin(&paf),
        expected(
            "join02",
            &["source=foo;idx=0;aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60"]
        )
    );
}

/// Run `badclip <args>` with stdin=/dev/null; assert it exits 2 and prints the
/// subcommand's usage (rather than blocking on stdin or erroring).
fn assert_shows_help(args: &[&str], usage: &str) {
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .args(args)
        .stdin(Stdio::null())
        .output()
        .expect("failed to run badclip");
    assert_eq!(output.status.code(), Some(2), "args={args:?}");
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(help.contains(usage), "expected {usage:?}, got: {help}");
}

#[test]
fn no_input_shows_help() {
    // Each subcommand with no file must print help and exit 2, not block on
    // stdin. stdin is /dev/null so a regression would EOF rather than hang.
    assert_shows_help(&["extract"], "Usage: badclip extract");
    assert_shows_help(&["geteseq"], "Usage: badclip geteseq");
    assert_shows_help(&["flteseq"], "Usage: badclip flteseq");
    // flteseq needs two inputs; one is still "no input" -> help, not an error.
    assert_shows_help(&["flteseq", "some.clip"], "Usage: badclip flteseq");
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
            &["source=foo;idx=0;aln_len=16881,2475;qlen=16848,15,2474;mapq=60,60"]
        )
    );
}

#[test]
fn paf_keeps_inversion() {
    // inv01.paf is a 3-hit inversion read (chr10 tp:A:P -, tp:A:I +, tp:A:P -).
    // Keeping the tp:A:I line yields two joins (the inverted middle segment),
    // matching what BAM produces from the SA tag.
    let expected = "\
chr10\t91443587\t<>\tchr10\t91446029\tm84039_230117_233243_s1/251662187/ccs\t60\t-\tsource=foo;idx=0;aln_len=2443,19097;qlen=2926,0,19078;mapq=60,60
chr10\t91443587\t><\tchr10\t91446029\tm84039_230117_233243_s1/251662187/ccs\t60\t-\tsource=foo;idx=1;aln_len=485,2443;qlen=485,0,21519;mapq=60,60
";
    assert_eq!(run_extract_file("inv01"), expected);
}

#[test]
fn source_tag_custom() {
    // -s overrides the default `source=foo` tag on every emitted record.
    let out = run_extract_file_args("join01", &["-s", "tumor"]);
    assert!(
        out.lines().all(|l| l
            .split('\t')
            .nth(8)
            .is_some_and(|info| info.split(';').any(|kv| kv == "source=tumor"))),
        "every record's INFO should carry source=tumor:\n{out}"
    );
    // The default is `source=foo`.
    assert!(run_extract_file("join01").contains("source=foo;"));
}

// --- CRAM input ---

// cram01.{bam,cram} are the SAME two synthetic reads (a right-soft-clip and a
// two-hit chimera) over the synthetic reference test/cram01.fa. CRAM stores the
// sequence as diffs against the reference, so `extract` must decode it via -r.

#[test]
fn cram_matches_bam() {
    // The autodetected CRAM path (with -r) must reproduce the BAM output exactly,
    // including the reference-reconstructed eseq bases.
    let bam = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg(test_dir().join("cram01.bam"))
        .output()
        .expect("failed to run badclip on BAM");
    assert!(bam.status.success());

    let cram = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg(test_dir().join("cram01.cram"))
        .arg("-r")
        .arg(test_dir().join("cram01.fa"))
        .output()
        .expect("failed to run badclip on CRAM");
    assert!(
        cram.status.success(),
        "badclip failed on CRAM: {}",
        String::from_utf8_lossy(&cram.stderr)
    );

    let cram_out = String::from_utf8(cram.stdout).unwrap();
    assert_eq!(cram_out, String::from_utf8(bam.stdout).unwrap());
    // cram01 has constant QUAL Q40, so `equal` is 40 and sits between elen and
    // eseq (the reference-decoded sequence carries usable qualities under CRAM).
    assert!(
        cram_out.contains(";equal=40;eseq="),
        "expected equal=40 before eseq:\n{cram_out}"
    );
}

#[test]
fn cram_requires_reference() {
    // CRAM without -r must fail with a message naming the reference / -r, not
    // emit a cryptic decode error.
    let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg(test_dir().join("cram01.cram"))
        .output()
        .expect("failed to run badclip");
    assert!(!out.status.success(), "CRAM without -r should fail");
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert!(
        stderr.contains("-r") && stderr.to_lowercase().contains("reference"),
        "expected a reference/-r error, got: {stderr}"
    );
}

#[test]
fn no_qual_omits_equal() {
    // A record may carry SEQ but no base qualities (QUAL = `*`). Then `eseq` is
    // still emitted but `equal` cannot be computed, so it is omitted. Fed as an
    // inline SAM on stdin (format is autodetected). 200M60S -> a right clip whose
    // 60 bp exceeds the default -c, so a clip breakend with eseq is emitted.
    let seq = "ACGT".repeat(65); // 260 bp = 200 aligned + 60 clipped
    let sam = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:400\nr1\t0\tchrA\t1\t60\t200M60S\t*\t0\t0\t{seq}\t*\n"
    );
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(sam.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(stdout.contains(";eseq="), "eseq should be emitted:\n{stdout}");
    assert!(
        !stdout.contains("equal="),
        "equal must be omitted when QUAL is absent:\n{stdout}"
    );
}

// --- BAM input ---

// bam01.bam holds three reads: a 3-hit chimera (.../234884627/ccs -> 1 left clip
// + 2 joins whose gaps exceed -e, so eseq is dropped but elen kept), a
// single-alignment read with a long clip (.../234884533/ccs), and a
// reverse-strand 2-hit translocation (.../240521013/ccs, join eseq present ->
// exercises reverse-complementing). Expected output is the committed golden file
// test/bam01.expected (eseq bytes verified against samtools).

#[test]
fn bam_matches() {
    let expected = std::fs::read_to_string(test_dir().join("bam01.expected")).unwrap();
    assert_eq!(run_extract_bam("bam01"), expected);
}

#[test]
fn bam_sorted_equals_unsorted() {
    // The SA-tag approach needs no grouping, so a coordinate-sorted BAM yields
    // the same breakend set (different line order).
    let mut unsorted: Vec<String> = run_extract_bam("bam01").lines().map(String::from).collect();
    let mut sorted: Vec<String> = run_extract_bam("bam01.srt")
        .lines()
        .map(String::from)
        .collect();
    unsorted.sort();
    sorted.sort();
    assert_eq!(unsorted, sorted);
}

#[test]
fn bam_paf_parity_single_hit() {
    // For a single-alignment read (no supplementary hits), BAM and PAF agree on
    // everything up to the BAM-only elen/eseq tags. bam01.paf is that read's PAF
    // line; strip the elen/eseq suffix from the BAM line before comparing.
    let paf_out = run_extract_file("bam01");
    let bam_line = run_extract_bam("bam01")
        .lines()
        .find(|l| l.contains("234884533"))
        .expect("single-hit read missing from BAM output")
        .to_string();
    let bam_trimmed = bam_line.split(";elen=").next().unwrap();
    assert_eq!(paf_out.trim_end(), bam_trimmed);
}

/// Parse the `elen=l,q,r` and optional `eseq=..` from a record's INFO column.
fn parse_elen_eseq(line: &str) -> Option<((i64, i64, i64), Option<String>)> {
    let info = line.split('\t').nth(8)?;
    let mut elen = None;
    let mut eseq = None;
    for kv in info.split(';') {
        if let Some(v) = kv.strip_prefix("elen=") {
            let n: Vec<i64> = v.split(',').map(|x| x.parse().unwrap()).collect();
            elen = Some((n[0], n[1], n[2]));
        } else if let Some(v) = kv.strip_prefix("eseq=") {
            eseq = Some(v.to_string());
        }
    }
    elen.map(|e| (e, eseq))
}

#[test]
fn bam_eseq_invariants() {
    // Every BAM line has elen; eseq is present iff the window (left+|qdist|+right)
    // is within the default -e (1000), and then its length equals that window.
    let out = run_extract_bam("bam01");
    let mut saw_present = false;
    let mut saw_dropped = false;
    for line in out.lines() {
        let ((l, q, r), eseq) = parse_elen_eseq(line).expect("BAM line missing elen");
        let window = l + q.abs() + r;
        match eseq {
            Some(s) => {
                assert!(window <= 1000, "eseq present but window {window} > 1000");
                assert_eq!(s.len() as i64, window, "eseq length != elen window");
                saw_present = true;
            }
            None => {
                assert!(window > 1000, "eseq dropped but window {window} <= 1000");
                saw_dropped = true;
            }
        }
    }
    assert!(
        saw_present && saw_dropped,
        "fixture should cover both cases"
    );
}

#[test]
fn bam_eseq_dropped_when_over_limit() {
    // With -e 0 every window exceeds the limit, so no eseq is emitted, but elen
    // is still present on every line.
    let out = run_extract_bam_args("bam01", &["-e", "0"]);
    for line in out.lines() {
        let (_, eseq) = parse_elen_eseq(line).expect("BAM line missing elen");
        assert!(eseq.is_none(), "eseq should be dropped with -e 0: {line}");
    }
}

// --- geteseq ---

/// Run `badclip geteseq -` feeding `bytes` on stdin.
fn run_geteseq_stdin(bytes: &[u8]) -> String {
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("geteseq")
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip");
    child.stdin.take().unwrap().write_all(bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "badclip geteseq failed");
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn geteseq_name_format_and_skips() {
    // First record has eseq -> FASTA name readName_idx_leftFlank_rightFlank
    // (idx and elen[0],elen[2]); second lacks eseq -> skipped.
    let input = "\
chrA\t100\t>>\tchrB\t200\tread1\t60\t+\tidx=2;aln_len=10,20;qlen=5,3,7;mapq=60,60;elen=5,3,7;eseq=ACGTACGTACGTAAA
chrC\t9\t<.\t.\t.\tread2\t30\t-\tidx=0;aln_len=8,0;qlen=8,0,2;mapq=30,0;elen=8,0,2
";
    assert_eq!(
        run_geteseq_stdin(input.as_bytes()),
        ">read1_2_5_7\nACGTACGTACGTAAA\n"
    );
}

#[test]
fn geteseq_matches_golden() {
    // geteseq over the extract golden output yields the committed FASTA.
    let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("geteseq")
        .arg(test_dir().join("bam01.expected"))
        .output()
        .expect("failed to run badclip geteseq");
    assert!(out.status.success());
    let expected = std::fs::read_to_string(test_dir().join("bam01.fa")).unwrap();
    assert_eq!(String::from_utf8(out.stdout).unwrap(), expected);
}

// --- flteseq ---

#[test]
fn flteseq_filters() {
    // flt01: elen=250,0,250 -> protected interval [200,300] with -l 50.
    // readA: no eseq -> dropped. readB: aln [0,500] contains [200,300] -> dropped.
    // readC: aln [250,500] does not contain -> kept. readD: no PAF line -> kept.
    // readE: two alns, one ([0,400]) contains -> dropped.
    let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("flteseq")
        .arg(test_dir().join("flt01.clip"))
        .arg(test_dir().join("flt01.rb3.paf"))
        .output()
        .expect("failed to run badclip flteseq");
    assert!(out.status.success());
    let expected = "\
chrD\t4\t<<\tchrE\t5\treadC\t60\t-\tsource=foo;idx=0;aln_len=1,1;qlen=1,0,1;mapq=60,60;elen=250,0,250;eseq=CCCC
chrF\t6\t>.\t.\t.\treadD\t60\t+\tsource=foo;idx=0;aln_len=1,0;qlen=1,0,1;mapq=60,0;elen=250,0,250;eseq=GGGG
";
    assert_eq!(String::from_utf8(out.stdout).unwrap(), expected);
}

#[test]
fn flteseq_source_relabel() {
    // With `-s novel`, all 5 input lines are printed; the survivors (readC, readD)
    // get their source rewritten to `novel`, and the dropped lines (readA no eseq,
    // readB and readE pangenome-explained) keep their original `source=foo`.
    let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("flteseq")
        .arg("-s")
        .arg("novel")
        .arg(test_dir().join("flt01.clip"))
        .arg(test_dir().join("flt01.rb3.paf"))
        .output()
        .expect("failed to run badclip flteseq");
    assert!(out.status.success());
    let expected = "\
chrA\t1\t>.\t.\t.\treadA\t60\t+\tsource=foo;idx=0;aln_len=100,0;qlen=100,0,60;mapq=60,0;elen=100,0,60
chrB\t2\t>>\tchrC\t3\treadB\t60\t+\tsource=foo;idx=0;aln_len=1,1;qlen=1,0,1;mapq=60,60;elen=250,0,250;eseq=AAAA
chrD\t4\t<<\tchrE\t5\treadC\t60\t-\tsource=novel;idx=0;aln_len=1,1;qlen=1,0,1;mapq=60,60;elen=250,0,250;eseq=CCCC
chrF\t6\t>.\t.\t.\treadD\t60\t+\tsource=novel;idx=0;aln_len=1,0;qlen=1,0,1;mapq=60,0;elen=250,0,250;eseq=GGGG
chrG\t7\t>>\tchrH\t8\treadE\t60\t+\tsource=foo;idx=0;aln_len=1,1;qlen=1,0,1;mapq=60,60;elen=250,0,250;eseq=TTTT
";
    assert_eq!(String::from_utf8(out.stdout).unwrap(), expected);
}
