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

/// Run `badclip extract [args...] -` feeding an inline SAM on stdin (the format
/// is autodetected) and return its stdout.
fn run_extract_sam(sam: &str, args: &[&str]) -> String {
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .args(args)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(sam.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(
        out.status.success(),
        "badclip failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn clip01_matches() {
    // Single hit (alen=14732); left trim (qs=141) then right trim (qlen-qe=292),
    // reported only with -T.
    assert_eq!(
        run_extract_file_args("clip01", &["-T"]),
        expected(
            "clip01",
            &[
                "source=foo;ctype=trim;idx=0;n_aln=1;aln_len=14732,0;qlen=14925,0,141;mapq=1,0",
                "source=foo;ctype=trim;idx=1;n_aln=1;aln_len=14732,0;qlen=14774,0,292;mapq=1,0",
            ]
        )
    );
    // Trims are off by default.
    assert_eq!(run_extract_file("clip01"), "");
}

#[test]
fn join01_matches() {
    // Flipped: chr22 (y1) first, chrY (y0) second.
    assert_eq!(
        run_extract_file("join01"),
        expected(
            "join01",
            &[
                "source=foo;ctype=split;idx=0;n_aln=2;aln_len=16881,2475;qlen=16848,15,2474;mapq=60,60"
            ]
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
            &[
                "source=foo;ctype=split;idx=0;n_aln=2;aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60"
            ]
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
            &[
                "source=foo;ctype=split;idx=0;n_aln=2;aln_len=29505,5668;qlen=29436,1,5661;mapq=60,60"
            ]
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
    // then a lone alignment, so only its left trim (qs=2489) is emitted (-T).
    assert_eq!(
        run_extract_file_args("join01", &["-T", "-a", "3000"]),
        "chr22\t22131975\t<.\t.\t.\tm84039_230117_233243_s1/257233489/ccs\t60\t+\
         \tsource=foo;ctype=trim;idx=0;n_aln=1;aln_len=16881,0;qlen=16848,0,2489;mapq=60,0\n"
    );
}

#[test]
fn stdin_matches_file() {
    let paf = std::fs::read(test_dir().join("join02.paf")).unwrap();
    assert_eq!(
        run_extract_stdin(&paf),
        expected(
            "join02",
            &[
                "source=foo;ctype=split;idx=0;n_aln=2;aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60"
            ]
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
    assert_shows_help(&["getfa"], "Usage: badclip getfa");
    assert_shows_help(&["fltbnd"], "Usage: badclip fltbnd");
    // fltbnd needs two inputs; one is still "no input" -> help, not an error.
    assert_shows_help(&["fltbnd", "some.clip"], "Usage: badclip fltbnd");
    // The pre-0.2 names still work as aliases.
    assert_shows_help(&["geteseq"], "Usage: badclip getfa");
    assert_shows_help(&["flteseq"], "Usage: badclip fltbnd");
}

#[test]
fn old_names_hidden_from_help() {
    let output = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("--help")
        .output()
        .expect("failed to run badclip");
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(help.contains("getfa") && help.contains("fltbnd"), "{help}");
    assert!(
        !help.contains("geteseq") && !help.contains("flteseq"),
        "{help}"
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
            &[
                "source=foo;ctype=split;idx=0;n_aln=2;aln_len=16881,2475;qlen=16848,15,2474;mapq=60,60"
            ]
        )
    );
}

#[test]
fn paf_keeps_inversion() {
    // inv01.paf is a 3-hit inversion read (chr10 tp:A:P -, tp:A:I +, tp:A:P -).
    // Keeping the tp:A:I line yields two joins (the inverted middle segment),
    // matching what BAM produces from the SA tag.
    let expected = "\
chr10\t91443587\t<>\tchr10\t91446029\tm84039_230117_233243_s1/251662187/ccs\t60\t-\tsource=foo;ctype=split;idx=0;n_aln=3;aln_len=2443,19097;qlen=2926,0,19078;mapq=60,60
chr10\t91443587\t><\tchr10\t91446029\tm84039_230117_233243_s1/251662187/ccs\t60\t-\tsource=foo;ctype=split;idx=1;n_aln=3;aln_len=485,2443;qlen=485,0,21519;mapq=60,60
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
        .arg("-T")
        .arg(test_dir().join("cram01.bam"))
        .output()
        .expect("failed to run badclip on BAM");
    assert!(bam.status.success());

    let cram = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("extract")
        .arg("-T")
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
    // 60 bp exceeds the default -c, so (with -T) a trim breakend with eseq is
    // emitted.
    let seq = "ACGT".repeat(65); // 260 bp = 200 aligned + 60 clipped
    let sam = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:400\nr1\t0\tchrA\t1\t60\t200M60S\t*\t0\t0\t{seq}\t*\n"
    );
    let stdout = run_extract_sam(&sam, &["-T"]);
    assert!(
        stdout.contains(";eseq="),
        "eseq should be emitted:\n{stdout}"
    );
    assert!(
        !stdout.contains("equal="),
        "equal must be omitted when QUAL is absent:\n{stdout}"
    );
}

#[test]
fn extract_min_mapq_is_a_col7_postfilter() {
    // -q is a post-filter on the emitted line's col-7 mapq; it does NOT drop
    // alignments upfront. rclip is a lone alignment (mapq 10 -> clip col7=10).
    // rjoin is a chimera whose primary is mapq 60 and SA supplementary mapq 5, so
    // the join's col7 = max(60,5) = 60. Fed as an inline SAM on stdin.
    let s160 = "ACGT".repeat(40);
    let q160 = "I".repeat(160);
    let s200 = "ACGT".repeat(50);
    let q200 = "I".repeat(200);
    let s100 = "ACGT".repeat(25);
    let q100 = "I".repeat(100);
    let sam = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:400\n@SQ\tSN:chrB\tLN:400\n\
         rclip\t0\tchrA\t1\t10\t100M60S\t*\t0\t0\t{s160}\t{q160}\n\
         rjoin\t0\tchrA\t1\t60\t100M100S\t*\t0\t0\t{s200}\t{q200}\tSA:Z:chrB,1,+,100S100M,5,0\n\
         rjoin\t2048\tchrB\t1\t5\t100H100M\t*\t0\t0\t{s100}\t{q100}\tSA:Z:chrA,1,+,100M100S,60,0\n"
    );
    // -T: the clip-shaped reads only emit trim breakends.
    let run = |args: &[&str]| -> String { run_extract_sam(&sam, &[&["-T"], args].concat()) };

    // Default -q 0: both the col7=10 clip and the col7=60 join are emitted.
    let all = run(&[]);
    assert!(
        all.contains("rclip") && all.contains("rjoin"),
        "default:\n{all}"
    );
    // -q 20: the clip (col7=10) is dropped; the join (col7=60) survives, and its
    // mapq tag still shows the mapq-5 hit -> the alignment was NOT dropped upfront.
    let filtered = run(&["-q", "20"]);
    assert!(
        !filtered.contains("rclip"),
        "clip col7=10 should drop:\n{filtered}"
    );
    assert!(
        filtered.contains("rjoin")
            && filtered.contains("n_aln=2;aln_len=100,100")
            && filtered.contains("mapq=60,5"),
        "join col7=60 should survive with both alignments intact:\n{filtered}"
    );
}

#[test]
fn extract_alt_contigs() {
    // --alt lists ALT contigs. A read whose primary lands on an ALT contig is
    // dropped; an ALT hit in the SA tag is filtered before clips/joins form.
    // `altprim` maps to chrAlt; `normal` maps to chrA with an SA to chrAlt.
    let s160 = "ACGT".repeat(40);
    let q160 = "I".repeat(160);
    let s200 = "ACGT".repeat(50);
    let q200 = "I".repeat(200);
    let sam = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:400\n@SQ\tSN:chrAlt\tLN:400\n\
         altprim\t0\tchrAlt\t1\t60\t100M60S\t*\t0\t0\t{s160}\t{q160}\n\
         normal\t0\tchrA\t1\t60\t100M100S\t*\t0\t0\t{s200}\t{q200}\tSA:Z:chrAlt,1,+,100S100M,60,0\n"
    );
    // A bwa-kit `.alt`-style file: `@` header lines are ignored, and the first
    // (tab-delimited) column of each record is the ALT contig name.
    let alt_path = format!("{}/alt_contigs.txt", env!("CARGO_TARGET_TMPDIR"));
    std::fs::write(
        &alt_path,
        "@HD\tVN:1.6\n@SQ\tSN:chrA\tLN:400\nchrAlt\t0\tchrA\t1\t60\t100M\t*\t0\t0\t*\t*\n",
    )
    .unwrap();

    // -T: the clip-shaped reads only emit trim breakends.
    let run = |args: &[&str]| -> String { run_extract_sam(&sam, &[&["-T"], args].concat()) };

    // Without --alt: altprim clips on chrAlt, normal joins chrA<->chrAlt.
    let plain = run(&[]);
    assert!(
        plain.contains("altprim") && plain.contains("chrAlt"),
        "plain:\n{plain}"
    );

    // With --alt: altprim (primary on ALT) is gone; normal's SA-to-ALT is dropped,
    // leaving only a chrA clip -> chrAlt appears nowhere.
    let filtered = run(&["--alt", &alt_path]);
    assert!(
        !filtered.contains("chrAlt") && !filtered.contains("altprim"),
        "ALT reads/hits should be excluded:\n{filtered}"
    );
    // `normal` survives as a plain clip on chrA (a `>.` clip, not a join).
    let cols: Vec<&str> = filtered.trim_end().split('\t').collect();
    assert_eq!(cols[0], "chrA");
    assert_eq!(cols[2], ">.");
    assert_eq!(cols[5], "normal");
}

#[test]
fn extract_min_equal_drops_low_quality() {
    // -Q (default 0 = keep all) drops breakends whose computed `equal` is below
    // the threshold, to shrink downstream input. rlow's window is all phred 5
    // (equal 5), rhigh's all phred 40 (equal 40); each 100M60S read emits one
    // right clip. Fed as an inline SAM on stdin.
    let seq = "ACGT".repeat(40); // 160 bp = 100 aligned + 60 clipped
    let qlow = "&".repeat(160); // ASCII 38 -> phred 5
    let qhigh = "I".repeat(160); // ASCII 73 -> phred 40
    let sam = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:400\n\
         rlow\t0\tchrA\t1\t60\t100M60S\t*\t0\t0\t{seq}\t{qlow}\n\
         rhigh\t0\tchrA\t1\t60\t100M60S\t*\t0\t0\t{seq}\t{qhigh}\n"
    );

    // -T: the clip-shaped reads only emit trim breakends.
    let run = |args: &[&str]| -> String { run_extract_sam(&sam, &[&["-T"], args].concat()) };

    // Default -Q 0: both breakends emitted (with their equal tags).
    let all = run(&[]);
    assert!(
        all.contains("rlow") && all.contains("rhigh"),
        "default:\n{all}"
    );
    // -Q 20: rlow (equal 5) dropped; rhigh (equal 40) kept.
    let filtered = run(&["-Q", "20"]);
    assert!(
        filtered.contains("rhigh") && !filtered.contains("rlow"),
        "rlow (equal=5) should be dropped by -Q 20:\n{filtered}"
    );
}

#[test]
fn paired_end_read_suffix() {
    // For paired-end reads the mates are disambiguated by a /1 (first segment,
    // FLAG 0x40) or /2 (last segment, FLAG 0x80) suffix on the read name; an
    // unpaired read (FLAG without 0x1) keeps its bare name. Each 100M60S read has
    // a 60 bp right clip (> default -c) so it emits one trim breakend (-T).
    let seq = "ACGT".repeat(40); // 160 bp = 100 aligned + 60 clipped
    let sam = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:400\n\
         mate\t65\tchrA\t1\t60\t100M60S\t=\t1\t0\t{seq}\t*\n\
         mate\t129\tchrA\t1\t60\t100M60S\t=\t1\t0\t{seq}\t*\n\
         solo\t0\tchrA\t1\t60\t100M60S\t*\t0\t0\t{seq}\t*\n"
    );
    let names: Vec<String> = run_extract_sam(&sam, &["-T"])
        .lines()
        .map(|l| l.split('\t').nth(5).unwrap().to_string())
        .collect();
    assert_eq!(names, vec!["mate/1", "mate/2", "solo"]);
}

// --- BAM input ---

// bam01.bam holds three reads: a 3-hit chimera (.../234884627/ccs -> 1 left trim
// + 2 splits whose gaps exceed -e, so eseq is dropped but elen kept), a
// single-alignment read with a long clip (.../234884533/ccs), and a
// reverse-strand 2-hit translocation (.../240521013/ccs, split eseq present ->
// exercises reverse-complementing). Expected output is the committed golden file
// test/bam01.expected, made with -T (eseq bytes verified against samtools).

#[test]
fn bam_matches() {
    let expected = std::fs::read_to_string(test_dir().join("bam01.expected")).unwrap();
    assert_eq!(run_extract_bam_args("bam01", &["-T"]), expected);
}

#[test]
fn trims_off_by_default() {
    // Without -T the trim lines are simply absent; every other line (and its
    // idx: a hidden trim keeps its number) is unchanged.
    let golden = std::fs::read_to_string(test_dir().join("bam01.expected")).unwrap();
    let expected: String = golden
        .lines()
        .filter(|l| !l.contains(";ctype=trim;"))
        .map(|l| format!("{l}\n"))
        .collect();
    let got = run_extract_bam("bam01");
    assert_eq!(got, expected);
    assert!(
        got.contains(";ctype=split;idx=1;"),
        "hidden trim keeps idx 0:\n{got}"
    );
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
    let paf_out = run_extract_file_args("bam01", &["-T"]);
    let bam_line = run_extract_bam_args("bam01", &["-T"])
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
    let out = run_extract_bam_args("bam01", &["-T"]);
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
    let out = run_extract_bam_args("bam01", &["-T", "-e", "0"]);
    for line in out.lines() {
        let (_, eseq) = parse_elen_eseq(line).expect("BAM line missing elen");
        assert!(eseq.is_none(), "eseq should be dropped with -e 0: {line}");
    }
}

// --- getfa ---

/// Run `badclip getfa -` feeding `bytes` on stdin.
fn run_getfa_stdin(bytes: &[u8]) -> String {
    let mut child = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("getfa")
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn badclip");
    child.stdin.take().unwrap().write_all(bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "badclip getfa failed");
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn getfa_name_format_and_skips() {
    // A record with eseq -> FASTA name readName_idx_leftFlank_rightFlank_code
    // (idx, elen[0], elen[2], and T/S/D/L/R for trim/split/del/ins); one without
    // eseq -> skipped. An insertion gives two records: L (eseq) then R (eseq2).
    let input = "\
chrA\t100\t>>\tchrB\t200\tread1\t60\t+\tctype=split;idx=2;aln_len=10,20;qlen=5,3,7;mapq=60,60;elen=5,3,7;eseq=ACGTACGTACGTAAA
chrC\t9\t<.\t.\t.\tread2\t30\t-\tctype=trim;idx=0;aln_len=8,0;qlen=8,0,2;mapq=30,0;elen=8,0,2
chrC\t9\t<.\t.\t.\tread3\t30\t-\tctype=trim;idx=1;elen=2,0,2;eseq=GGTT
chrD\t50\t>>\tchrD\t1550\tread4\t60\t+\tctype=del;idx=300;elen=2,0,2;eseq=CCAA
chrE\t70\t>>\tchrE\t74\tread5\t60\t-\tctype=ins;idx=400;inslen=1000;elen=2,0,1;eseq=AAC;elen2=1,0,2;eseq2=GTT
";
    assert_eq!(
        run_getfa_stdin(input.as_bytes()),
        ">read1_2_5_7_S\nACGTACGTACGTAAA\n>read3_1_2_2_T\nGGTT\n>read4_300_2_2_D\nCCAA\n\
         >read5_400_2_1_L\nAAC\n>read5_400_1_2_R\nGTT\n"
    );
}

#[test]
fn getfa_infers_missing_ctype() {
    // Pre-0.2 extract output has no ctype: a mate-less line is a trim, else a split.
    let input = "\
chrA\t1\t>.\t.\t.\tr1\t60\t+\tidx=0;elen=1,0,1;eseq=AC
chrA\t1\t>>\tchrB\t2\tr2\t60\t+\tidx=0;elen=1,0,1;eseq=GT
";
    assert_eq!(
        run_getfa_stdin(input.as_bytes()),
        ">r1_0_1_1_T\nAC\n>r2_0_1_1_S\nGT\n"
    );
}

#[test]
fn getfa_min_equal_filters() {
    // -Q (default 20) drops records whose `equal` is below the threshold; records
    // without an `equal` tag are kept.
    let input = "\
chrA\t1\t>.\t.\t.\trlow\t60\t+\tidx=0;elen=5,0,7;equal=10;eseq=AAAA
chrB\t2\t>.\t.\t.\trhigh\t60\t+\tidx=0;elen=5,0,7;equal=30;eseq=CCCC
chrC\t3\t>.\t.\t.\trnoq\t60\t+\tidx=0;elen=5,0,7;eseq=GGGG
chrD\t4\t>>\tchrD\t4\trins\t60\t+\tctype=ins;idx=9;elen=1,0,1;equal=30;eseq=AC;elen2=1,0,1;equal2=10;eseq2=GT
";
    // Default -Q 20: rlow (equal=10) dropped; rhigh (equal=30) and rnoq (no
    // equal) kept; rins is dropped whole, as its lower window quality is 10.
    assert_eq!(
        run_getfa_stdin(input.as_bytes()),
        ">rhigh_0_5_7_T\nCCCC\n>rnoq_0_5_7_T\nGGGG\n"
    );
}

#[test]
fn getfa_matches_golden() {
    // getfa over the extract golden output yields the committed FASTA; the
    // `geteseq` alias behaves identically.
    let expected = std::fs::read_to_string(test_dir().join("bam01.fa")).unwrap();
    for cmd in ["getfa", "geteseq"] {
        let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
            .arg(cmd)
            .arg(test_dir().join("bam01.expected"))
            .output()
            .expect("failed to run badclip getfa");
        assert!(out.status.success());
        assert_eq!(String::from_utf8(out.stdout).unwrap(), expected, "{cmd}");
    }
}

// --- fltbnd ---

#[test]
fn fltbnd_filters() {
    // flt01: elen=250,0,250 -> protected interval [200,300] with -l 50.
    // readA: no eseq -> dropped. readB: aln [0,500] contains [200,300] -> dropped.
    // readC: aln [250,500] does not contain -> kept. readD: no PAF line -> kept.
    // readE: two alns, one ([0,400]) contains -> dropped.
    let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("fltbnd")
        .arg(test_dir().join("flt01.clip"))
        .arg(test_dir().join("flt01.rb3.paf"))
        .output()
        .expect("failed to run badclip fltbnd");
    assert!(out.status.success());
    let expected = "\
chrD\t4\t<<\tchrE\t5\treadC\t60\t-\tsource=foo;idx=0;aln_len=1,1;qlen=1,0,1;mapq=60,60;elen=250,0,250;eseq=CCCC
chrF\t6\t>.\t.\t.\treadD\t60\t+\tsource=foo;idx=0;aln_len=1,0;qlen=1,0,1;mapq=60,0;elen=250,0,250;eseq=GGGG
";
    assert_eq!(String::from_utf8(out.stdout).unwrap(), expected);
}

#[test]
fn fltbnd_source_relabel() {
    // With `-s novel`, all 5 input lines are printed; the survivors (readC, readD)
    // get their source rewritten to `novel`, and the dropped lines (readA no eseq,
    // readB and readE pangenome-explained) keep their original `source=foo`.
    let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
        .arg("fltbnd")
        .arg("-s")
        .arg("novel")
        .arg(test_dir().join("flt01.clip"))
        .arg(test_dir().join("flt01.rb3.paf"))
        .output()
        .expect("failed to run badclip fltbnd");
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

#[test]
fn fltbnd_min_equal_filters() {
    // -Q drops lines whose `equal` is below the threshold. Two survivors (empty
    // PAF -> nothing protected); rlow has equal=5, rkeep equal=30.
    let dir = env!("CARGO_TARGET_TMPDIR");
    let clip = format!("{dir}/fltbnd_q.clip");
    let paf = format!("{dir}/fltbnd_q_empty.paf");
    std::fs::write(&paf, "").unwrap(); // empty PAF -> no line is "protected"
    std::fs::write(
        &clip,
        "chrA\t1\t>.\t.\t.\trkeep\t60\t+\tsource=foo;idx=0;elen=250,0,250;equal=30;eseq=AAAA\n\
         chrB\t2\t>.\t.\t.\trlow\t60\t+\tsource=foo;idx=0;elen=250,0,250;equal=5;eseq=CCCC\n",
    )
    .unwrap();

    let run = |args: &[&str]| -> String {
        let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
            .arg("fltbnd")
            .args(args)
            .arg(&clip)
            .arg(&paf)
            .output()
            .expect("failed to run badclip fltbnd");
        assert!(out.status.success());
        String::from_utf8(out.stdout).unwrap()
    };

    // Default -Q 20: both are survivors, but rlow (equal=5) is filtered out.
    let got = run(&[]);
    assert!(
        got.contains("rkeep") && !got.contains("rlow"),
        "rlow (equal=5) should be dropped by default -Q 20:\n{got}"
    );
    // -Q 0 keeps rlow too.
    assert!(
        run(&["-Q", "0"]).contains("rlow"),
        "rlow should survive -Q 0"
    );
}

#[test]
fn fltbnd_insertion_either_junction() {
    // An insertion line has two eseq windows, named _L and _R by getfa; it is
    // dropped if EITHER junction is spanned by a pangenome alignment. With
    // elen=250,0,250 and -l 50 each junction interval is [200,300].
    let dir = env!("CARGO_TARGET_TMPDIR");
    let clip = format!("{dir}/fltbnd_ins.clip");
    let paf = format!("{dir}/fltbnd_ins.paf");
    let w = "elen=250,0,250;eseq=AAAA;elen2=250,0,250;eseq2=CCCC";
    std::fs::write(
        &clip,
        format!(
            "chrA\t1\t>>\tchrA\t5\tinsL\t60\t+\tsource=foo;ctype=ins;idx=7;{w}\n\
             chrA\t1\t>>\tchrA\t5\tinsR\t60\t+\tsource=foo;ctype=ins;idx=7;{w}\n\
             chrA\t1\t>>\tchrA\t5\tinsN\t60\t+\tsource=foo;ctype=ins;idx=7;{w}\n"
        ),
    )
    .unwrap();
    // insL: its L junction is spanned. insR: only its R junction is spanned.
    // insN: both alignments stop short of the junction -> novel.
    std::fs::write(
        &paf,
        "insL_7_250_250_L\t500\t0\t500\t+\tPG\t1000\t0\t500\n\
         insL_7_250_250_R\t500\t250\t500\t+\tPG\t1000\t0\t250\n\
         insR_7_250_250_R\t500\t0\t500\t+\tPG\t1000\t0\t500\n\
         insN_7_250_250_L\t500\t250\t500\t+\tPG\t1000\t0\t250\n\
         insN_7_250_250_R\t500\t0\t250\t+\tPG\t1000\t0\t250\n",
    )
    .unwrap();
    let run = |cmd: &str, args: &[&str]| -> String {
        let out = Command::new(env!("CARGO_BIN_EXE_badclip"))
            .arg(cmd)
            .args(args)
            .arg(&clip)
            .arg(&paf)
            .output()
            .expect("failed to run badclip fltbnd");
        assert!(out.status.success());
        String::from_utf8(out.stdout).unwrap()
    };
    let got = run("fltbnd", &[]);
    assert_eq!(got.lines().count(), 1, "{got}");
    assert!(got.contains("\tinsN\t"), "only insN is novel:\n{got}");
    // -s prints everything, relabelling only insN; the flteseq alias agrees.
    let all = run("flteseq", &["-s", "novel"]);
    let novel: Vec<&str> = all.lines().filter(|l| l.contains("source=novel")).collect();
    assert_eq!(all.lines().count(), 3, "{all}");
    assert_eq!(novel.len(), 1, "{all}");
    assert!(novel[0].contains("\tinsN\t"), "{all}");
}

// --- CIGAR indels (del/ins) ---

/// Deterministic pseudo-random bases (an LCG), for synthetic reads.
fn bases(seed: u64, n: usize) -> String {
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 62) as usize] as char
        })
        .collect()
}

/// A base different from `c`.
fn other(c: char) -> char {
    if c == 'A' { 'C' } else { 'A' }
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            c => c,
        })
        .collect()
}

/// Synthetic reads exercising del/ins extraction (`-f 10` keeps eseq short),
/// returned as (SAM text, expected default output lines).
fn gap_sam() -> (String, Vec<String>) {
    let q = |n: usize| "I".repeat(n); // phred 40
    // g1: forward deletion. g2: reverse deletion (alignment-orientation q=80,
    // read-forward 200-80=120).
    let s1 = bases(1, 200);
    let s2 = bases(2, 200);
    let r2 = revcomp(&s2);
    // g3: forward insertion whose first 5 bases repeat the next read bases (a
    // TSD-like right shift of 5) and whose last 3 repeat the preceding bases (a
    // left shift of 3); the second window has low quality (phred 5).
    let (x3, y3) = (bases(3, 100), bases(4, 100));
    let mut i3: Vec<char> = bases(5, 1200).chars().collect();
    let (xc, yc): (Vec<char>, Vec<char>) = (x3.chars().collect(), y3.chars().collect());
    i3[..5].copy_from_slice(&yc[..5]);
    i3[5] = other(yc[5]);
    i3[1197..].copy_from_slice(&xc[97..]);
    i3[1196] = other(xc[96]);
    let s3 = format!("{x3}{}{y3}", i3.iter().collect::<String>());
    let q3 = format!("{}{}{}", q(1290), "&".repeat(20), q(90));
    // g4: reverse insertion (q=80, read-forward qa=1400-80-1200=120); on the
    // stored SEQ it shifts right by 2 and not left.
    let (x4, y4) = (bases(6, 80), bases(7, 120));
    let mut i4: Vec<char> = bases(8, 1200).chars().collect();
    let (xc, yc): (Vec<char>, Vec<char>) = (x4.chars().collect(), y4.chars().collect());
    i4[..2].copy_from_slice(&yc[..2]);
    i4[2] = other(yc[2]);
    i4[1199] = other(xc[79]);
    let s4 = format!("{x4}{}{y4}", i4.iter().collect::<String>());
    let r4 = revcomp(&s4);
    // g5: a chimera; its hard-clipped supplementary record (read [300,500))
    // carries a deletion 5 bp into its SEQ, so that eseq is clamped to 5+10.
    let s5 = bases(9, 500);
    let sam = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:10000\n@SQ\tSN:chrB\tLN:10000\n\
         g1\t0\tchrA\t1\t60\t100M1500D100M\t*\t0\t0\t{s1}\t{}\n\
         g2\t16\tchrA\t1\t60\t80M1500D120M\t*\t0\t0\t{s2}\t{}\n\
         g3\t0\tchrA\t1\t60\t100M1200I100M\t*\t0\t0\t{s3}\t{q3}\n\
         g4\t16\tchrA\t1\t60\t80M1200I120M\t*\t0\t0\t{s4}\t{}\n\
         g5\t0\tchrA\t2001\t60\t300M200S\t*\t0\t0\t{s5}\t{}\tSA:Z:chrB,1,+,300S200M1500D,60,0;\n\
         g5\t2048\tchrB\t1\t60\t300H5M1500D195M\t*\t0\t0\t{}\t{}\tSA:Z:chrA,2001,+,300M200S,60,0;\n",
        q(200),
        q(200),
        q(1400),
        q(500),
        &s5[300..],
        q(200),
    );
    let lines = vec![
        format!(
            "chrA\t100\t>>\tchrA\t1600\tg1\t60\t+\tsource=foo;ctype=del;idx=100;n_aln=1;\
             aln_len=1700,1700;qlen=100,0,100;mapq=60,60;elen=10,0,10;equal=40;eseq={}",
            &s1[90..110]
        ),
        format!(
            "chrA\t80\t>>\tchrA\t1580\tg2\t60\t-\tsource=foo;ctype=del;idx=120;n_aln=1;\
             aln_len=1700,1700;qlen=80,0,120;mapq=60,60;elen=10,0,10;equal=40;eseq={}",
            &r2[110..130]
        ),
        format!(
            "chrA\t97\t>>\tchrA\t105\tg3\t60\t+\tsource=foo;ctype=ins;idx=100;n_aln=1;\
             aln_len=1400,1400;qlen=100,1200,100;mapq=60,60;inslen=1200;\
             elen=10,0,10;equal=40;eseq={};elen2=10,0,10;equal2=5;eseq2={}",
            &s3[90..110],
            &s3[1290..1310]
        ),
        format!(
            "chrA\t80\t>>\tchrA\t82\tg4\t60\t-\tsource=foo;ctype=ins;idx=120;n_aln=1;\
             aln_len=1400,1400;qlen=80,1200,120;mapq=60,60;inslen=1200;\
             elen=10,0,10;equal=40;eseq={};elen2=10,0,10;equal2=40;eseq2={}",
            &r4[110..130],
            &r4[1310..1330]
        ),
        format!(
            "chrA\t2300\t>>\tchrB\t0\tg5\t60\t+\tsource=foo;ctype=split;idx=0;n_aln=2;\
             aln_len=300,1700;qlen=300,0,200;mapq=60,60;elen=10,0,10;equal=40;eseq={}",
            &s5[290..310]
        ),
        format!(
            "chrB\t5\t>>\tchrB\t1505\tg5\t60\t+\tsource=foo;ctype=del;idx=305;n_aln=2;\
             aln_len=1700,1700;qlen=305,0,195;mapq=60,60;elen=5,0,10;equal=40;eseq={}",
            &s5[300..315]
        ),
    ];
    (sam, lines)
}

/// Keep the lines of `expected` whose index is in `keep`, newline-terminated.
fn pick(expected: &[String], keep: &[usize]) -> String {
    keep.iter().map(|&i| format!("{}\n", expected[i])).collect()
}

#[test]
fn gap_del_ins_lines() {
    // Exact del/ins lines: positions (insertions normalized to their leftmost/
    // rightmost placement), strand, qlen, idx (read-forward offset), inslen,
    // and eseq/eseq2 windows (read-forward, clamped to a hard-clipped SEQ).
    let (sam, expected) = gap_sam();
    assert_eq!(
        run_extract_sam(&sam, &["-f", "10"]),
        pick(&expected, &[0, 1, 2, 3, 4, 5])
    );
}

#[test]
fn gap_min_gap() {
    // -g sets the minimum indel length (k suffix ok); -g 0 disables del/ins.
    let (sam, expected) = gap_sam();
    let run = |g: &str| run_extract_sam(&sam, &["-f", "10", "-g", g]);
    assert_eq!(
        run("1.5k"),
        pick(&expected, &[0, 1, 4, 5]),
        "-g 1.5k drops the 1200I"
    );
    assert_eq!(run("1500"), run("1.5k"));
    assert_eq!(
        run("2000"),
        pick(&expected, &[4]),
        "-g 2000 drops every indel"
    );
    assert_eq!(run("0"), pick(&expected, &[4]), "-g 0 disables del/ins");
}

#[test]
fn gap_min_equal_uses_both_windows() {
    // g3's second insertion window is phred 5: -Q 20 drops the whole line.
    let (sam, expected) = gap_sam();
    assert_eq!(
        run_extract_sam(&sam, &["-f", "10", "-Q", "20"]),
        pick(&expected, &[0, 1, 3, 4, 5])
    );
}

#[test]
fn gap_from_paf_cigar() {
    // PAF input takes del/ins from cg:Z. Without a read sequence there is no
    // eseq, and an insertion stays at its CIGAR placement (t >> t).
    let paf = "\
g1\t200\t0\t200\t+\tchrA\t10000\t0\t1700\t200\t1700\t60\ttp:A:P\tcg:Z:100M1500D100M
g2\t200\t0\t200\t-\tchrA\t10000\t0\t1700\t200\t1700\t60\ttp:A:P\tcg:Z:80M1500D120M
g4\t1400\t0\t1400\t-\tchrA\t10000\t0\t200\t200\t1400\t60\ttp:A:P\tcg:Z:80M1200I120M
";
    let (_, expected) = gap_sam();
    let strip = |l: &str| l.split(";elen=").next().unwrap().to_string();
    let want = format!(
        "{}\n{}\n{}\n",
        strip(&expected[0]),
        strip(&expected[1]),
        strip(&expected[3]).replace("chrA\t80\t>>\tchrA\t82", "chrA\t80\t>>\tchrA\t80")
    );
    assert_eq!(run_extract_stdin(paf.as_bytes()), want);
}

// gap01.bam holds every record of four HG002 HiFi reads with >=1 kb CIGAR
// indels: .../240521100/ccs (reverse primary insertion), .../240587166/ccs
// (reverse primary deletion), .../234886332/ccs (a 4-hit chimera whose
// hard-clipped reverse supplementary records carry an insertion and a
// deletion) and .../234886162/ccs (forward primary deletion). gap01.paf is the
// same reads' minimap2 PAF (cg:Z kept, cs:Z dropped). The golden output's
// positions, insertion shifts, qlen and eseq/equal values were checked
// independently against samtools view and the reads' FASTQ.

#[test]
fn gap01_matches() {
    let expected = std::fs::read_to_string(test_dir().join("gap01.expected")).unwrap();
    assert_eq!(run_extract_bam("gap01"), expected);
}

#[test]
fn gap01_bam_paf_parity() {
    // For the same alignments, PAF and BAM agree up to the read-sequence tags,
    // the aln_len of SA-derived split hits, and insertion normalization: the
    // PAF insertion sits at its CIGAR offset t, the BAM one spans [L, R] ∋ t.
    let norm = |out: &str| -> Vec<(String, Option<(i64, i64)>)> {
        let mut v: Vec<(String, Option<(i64, i64)>)> = out
            .lines()
            .map(|l| {
                let l = l.split(";elen=").next().unwrap();
                let mut f: Vec<String> = l.split('\t').map(String::from).collect();
                let mut span = None;
                if f[8].contains("ctype=ins") {
                    span = Some((f[1].parse().unwrap(), f[4].parse().unwrap()));
                    f[1] = ".".into();
                    f[4] = ".".into();
                } else if f[8].contains("ctype=split") {
                    let info: Vec<&str> = f[8]
                        .split(';')
                        .filter(|kv| !kv.starts_with("aln_len="))
                        .collect();
                    f[8] = info.join(";");
                }
                (f.join("\t"), span)
            })
            .collect();
        v.sort();
        v
    };
    let bam = norm(&run_extract_bam("gap01"));
    let paf = norm(&run_extract_file("gap01"));
    assert_eq!(bam.len(), paf.len());
    let mut n_ins = 0;
    for ((bl, bspan), (pl, pspan)) in bam.iter().zip(&paf) {
        assert_eq!(bl, pl);
        if let (Some((l, r)), Some((t, t2))) = (bspan, pspan) {
            assert_eq!(t, t2, "PAF insertion is not shifted");
            assert!(l <= t && t <= r, "{l} <= {t} <= {r}");
            n_ins += 1;
        }
    }
    assert_eq!(n_ins, 2, "fixture has two insertions");
}
