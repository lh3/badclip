# CLAUDE.md

Guidance for working in this repository.

## What this is

`badclip` is a Rust CLI that extracts breakends and structural-variant signals
from long-read alignments. It is a from-scratch reimplementation inspired by
`test/minisv.js` — a **reference, not ground truth** (it may be buggy, and we
implement only a simplified subset). The user's written spec and the `test/*`
fixtures are authoritative.

## Build & test

```sh
cargo build            # debug binary at target/debug/badclip
cargo test             # runs tests/extract.rs against the fixtures
cargo run -- extract test/bam01.bam          # BAM (default)
cargo run -- extract --paf test/join02.paf   # PAF
```

## Layout

- `src/main.rs`   — clap CLI, subcommand dispatch. New subcommands slot in here.
- `src/io.rs`     — `open_reader`: stdin/file + transparent gzip (magic-byte peek); PAF only.
- `src/paf.rs`    — `Hit` struct, `Strand`, `parse_paf_line` (keeps only `tp:A:P`).
- `src/aln.rs`    — SAM/BAM/CRAM input via noodles-util (format autodetected): build `Hit`s from a primary record + its `SA` tag. CRAM needs `-r`.
- `src/extract.rs`— dispatch (alignment file vs `--paf`), grouping, sorting, clip/join emission; shared `emit_read`/`passes_filter`.
- `src/geteseq.rs`— `geteseq` subcommand: `extract` output → FASTA of `eseq` records.
- `src/flteseq.rs`— `flteseq` subcommand: filter breakends by pangenome eseq alignment.
- `src/merge.rs`  — `merge` subcommand: cluster per-read breakends into consensus SV calls.
- `src/fltreg.rs` — `fltreg` subcommand: drop `extract`/`merge` breakends that fall in BED regions.
- `tests/extract.rs`, `tests/merge.rs`, `tests/fltreg.rs` — end-to-end tests driving the compiled binary.
- `test/`         — `*.paf` inputs, `*.msv` (minisv) references, `bam01*`/`cram01*`/`flt01*`/`merge01*` fixtures/goldens, `fltreg01.bed`, `minisv.js`.

## `geteseq`

Filters `extract` output (TAB-delimited; `-` = stdin, gzip auto-detected; no
input → print help) to FASTA. For each record with an `eseq` tag it emits
`>readName_idx_L_R` (from
the `idx` tag and `elen`'s `leftFlank`/`rightFlank`) followed by the `eseq` bases;
records lacking `eseq` are skipped. `-Q` (default 20) additionally drops records
whose `equal` quality is below the threshold; records without an `equal` tag are
kept. Golden: `test/bam01.fa`.

## `flteseq`

Filters `extract` output against a ropebwt3 `sw` PAF (both gzip-ok) of the
`geteseq` FASTA aligned to a pangenome. Drops a line if it has no `eseq`, or if
some PAF alignment's query interval `[qs,qe]` (cols 3,4) contains the junction
interval `[max(0,e0-l), min(e0+|e1|+l, e0+|e1|+e2)]` on the eseq (`elen=e0,e1,e2`,
`-l` margin, default 50) — the breakend is "protected"/not novel. Survivors are
printed verbatim. The PAF qname is `geteseq`'s name and its lines are in the same
order as the `eseq` lines (a subset), so both files stream with a one-line PAF
lookahead — no in-memory load. Either input may be `-` (stdin); missing input →
help (via `main.rs::print_subcommand_help`, shared with `extract`/`geteseq`).
Fixtures: `test/flt01.clip`, `test/flt01.rb3.paf`.

`-Q` (default 20) additionally drops any line whose `equal` quality is below the
threshold (in both modes; lines without an `equal` tag are kept). The check runs
**after** the PAF lookahead consumes that line's entries, so the stream stays in
sync regardless of how the PAF was built (`flteseq.rs::line_equal`).

`-s STR` flips the mode: print **all** input lines and rewrite the `source=` INFO
tag to `STR` on the survivors (kept/novel), leaving dropped lines (no `eseq`, or
pangenome-explained) verbatim — tagging novel calls distinctly for a downstream
`merge` (`flteseq.rs::relabel_source`). Without `-s`, only survivors are printed
verbatim.

## `merge`

Collapses per-read `extract` breakends into consensus SV calls — a simplified
`gc_cmd_merge` from `minisv.js`. badclip has only the breakend type (no
`SVTYPE`/`SVLEN`), no per-sample `source`, and no centromere/RT annotation, so
minisv's SVTYPE/SVLEN checks, the `-d` filter, the centromere filter (`-e`), and
the RT branch (`-r`/`-R`) are all dropped. Unlike minisv (which needs an upstream
`sort -k1,1 -k2,2n`), `merge` **loads all records into memory and sorts them
itself** by `(ctg, pos)` — Rust `str` order == `LC_ALL=C sort`; the sort is
stable so the representative pick is deterministic. Input `-`=stdin, gzip
auto-detected; no input → help.

Algorithm (`src/merge.rs`): sweep a window of active clusters. Each record joins
the active cluster whose members it most often matches (`same_sv`), else starts a
new one; a cluster is flushed (and emitted if it passes the filters) once the
sweep moves past its `pos_max + -w`. `same_sv` requires: same `ctg`; compatible
`ori` (equal, or the inversion pair `><`↔`<>`); `|pos−pos| ≤ -w`; same `ctg2`;
and `|pos2−pos2| ≤ -w` **skipped when either `pos2` is absent** (a clip's `.`
mate — reproduces minisv's `NaN` comparison, so clips cluster on `ctg2`/`pos`
alone). Per-cluster comparisons are capped at `-C` members (a deterministic cap
replacing minisv's reservoir sampling — deep-cluster counts are capped, not
scaled), and the active-cluster list is bounded by `-A`.

Flags mirror minisv's kept subset: `-c` min read count (3), `-s` min count on
each strand (1), `-w` window bp (100), `-A` max active clusters (100), `-M` max
reads compared per cluster (500). `-Q` (default 20) drops input breakends whose
`equal` quality is below the threshold before clustering (`parse_rec`); breakends
without an `equal` tag are kept. `-q` (default 0) drops input breakends whose
col-7 mapq (the join's `max`, or a clip's mapq) is below the threshold. `-p`
(default 10) drops input **join** breakends whose *smaller* per-hit mapq (the
`min` of the `mapq=` INFO pair) is below the threshold; clips are exempt (their
`mapq2` is 0, so their min would always be 0), so `-p` never affects a clip
line — `-q`/`-p` are the max/min poles of the same per-hit mapq pair
(`parse_rec`).

A cluster is emitted only if `count ≥ -c` and each strand has `≥ -s` reads. The
representative (`members[len/2]`) supplies the coordinate/`ori`/`strand` fields.
Output is 9 TAB columns, `extract`-shaped:

```
ctg1  pos1  ori  ctg2  pos2  .  count  strand  INFO
```

(`.` in the qname slot, cluster size in the mapq slot; `pos2`=`.` for a clip.)
`INFO` is merge-derived only (the representative's `idx`/`aln_len`/`qlen`/`mapq`/
`elen`/`eseq` are **not** carried through): `avg_mapq=q1,q2;count=<src:f,r|...>`
(`q1` = mean per-hit mapq on the `ctg1:pos1` side, `q2` on the `ctg2:pos2` side
from each member's `mapq=` INFO pair — already output-endpoint-aligned, so
inversion `><`/`<>` members average correctly; `q2=0` for a clip cluster) then,
when any member is `><`/`<>`, `;count_fr=A;count_rf=B` and `;foldback` if `A*B==0
&& ctg1==ctg2`, and finally `;reads=<src:name,...|...>`. `count=` lists each
`source=` observed in the cluster (from the reads' `source=` tag) with its
forward,reverse counts, `|`-joined and **alphabetical** by source; a source is
included only if present (never `x:0,0`), e.g. `count=foo:2,3|retain:4,5` or
single-source `count=foo:2,2`. `reads=` is stratified the same way — the
supporting read names per source, `|`-joined and alphabetical by source (names in
member order within a source), e.g. `reads=foo:r1,r2|retain:r3,r4` — so it lists
the same sources in the same order as `count=`. The `-s`/`-c` filters still use
cluster-wide strand totals.
**Caveat:** the col-8 strand is
`+`/`-` from join canonicalization (read orientation for clips), not strictly
read strand, so `-s` is an approximate two-sided-support heuristic. Fixtures:
`test/merge01.clip` (unsorted input), `test/merge01.expected` (default-threshold
golden).

### `-m` sample-merge

`-m` re-runs the **same clustering algorithm** on the (filtered) combination of
per-sample `merge` outputs to produce population-level, sample-merged calls.
`merge` and `extract` output share the 9 columns but differ in INFO, so `-m`
swaps only the input parser (`parse_rec_merge`) and output emitter
(`write_sv_merge`); `same_sv`/sweep are reused verbatim.

- **Input parsing** (`-m`): mapq from `avg_mapq=q1,q2` (col 7 is a count here,
  not mapq); per-strand read counts from `count=` summed over its
  `|`-separated `[src:]f,r` entries (source-prefix optional via `sum_count`, so
  `-m` output re-feeds cleanly); `count_fr`/`count_rf` read straight from those
  tags. `-Q` is a no-op (no `equal` tag). `-q`/`-p` apply to the
  `avg_mapq`-derived `max`/`min` (clips exempt from `-p`).
- **Per-line filters** (`-m` only, applied in `parse_rec_merge` before
  clustering): `-C` (default 3) drops input lines whose **total** count
  (Σ fwd+rev) is below the threshold; `-S` (default 1) drops lines whose
  per-strand count is below it on either strand.
- **Counting**: the `-c`/`-s` filters are **read-based** — each input line
  contributes its `count=` totals, so `-c` thresholds Σ(fwd+rev) and `-s` the
  per-strand read totals. But the **col-7 count slot is the number of distinct
  sources** (samples) in the cluster — the union of the `src:` labels across the
  members' `count=` tags (fallback: member count when the input carries no
  labels, e.g. a chained `-m` pass over `count=F,R`). `count_fr`/`count_rf` are
  summed from the input tags (read-weighted), `foldback` when
  `count_fr*count_rf==0 && ctg1==ctg2`. `avg_mapq=q1,q2` is the unweighted mean
  of members' per-end avg mapq.
- **Output INFO** (`-m`): `avg_mapq=q1,q2;count=F,R` (aggregate read `F,R`, **no**
  per-source breakdown) then, for inversions, `;count_fr=A;count_rf=B[;foldback]`.
  **No `reads=` tag** (too many at 3202-sample scale). The output is itself
  `-m`-parseable, so merges chain (col-7 then falls back to member count, since
  `count=F,R` drops the source labels).

The `-C`/`-S`/`-Q`/`-q`/`-p` per-line filters and `-c`/`-s` cluster filters are
all no-ops-or-unchanged in the default (non-`-m`) extract path.

## `fltreg`

Drops `extract`/`merge` output lines whose breakends land in BED regions.
`fltreg <in> <bed>` (both positional; `-` = stdin for `<in>`, gzip auto on both;
either missing → help). Since both output formats share the 9 columns, one filter
serves both: a line is dropped if **either** endpoint — `(ctg1,pos1)` (cols 0,1)
or, when present, `(ctg2,pos2)` (cols 3,4) — falls in a region; a clip
(`ctg2="."`) tests only endpoint 1. Survivors are printed verbatim; lines whose
endpoints can't be parsed are kept.

Positions are raw 0-based offsets and BED is half-open, so offset `p` is inside
`[start,end)` iff `start <= p < end`. BED loading (`Regions::load`) skips blank /
`#` / `track` / `browser` lines, takes cols 0/1/2 as `chrom start end`, and per
contig **sorts by start and merges overlapping/adjacent intervals** into a
disjoint list; `Regions::contains` then does one `partition_point` binary search
(the rightmost interval with `start <= p` is the only one that can contain `p`).
No new dependency. Fixture: `test/fltreg01.bed`.

## Interval / offset notation

An interval is a pair of **offsets** `|st,en|`, where an offset sits *between*
two adjacent bases, not on a base (0-based, half-open). In `"A|CGT|AGC"`,
`|1,4|` is `"CGT"`. A hit matches a target interval `tc:|ts,te|` to a query
interval `|qs,qe|` on a strand. Emitted positions are raw 0-based offsets — no
±1 adjustment.

## `extract` — behaviour

Input is an **alignment file by default** (SAM/BAM/CRAM, autodetected), or **PAF
with `--paf`**. `-` reads stdin (no input at all → print help, don't block on
stdin). Both paths build a `Vec<Hit>` per read and feed the shared `emit_read`.

- **PAF** (`src/extract.rs::run_paf`): assumed grouped by read name; streamed and
  buffered per read, keeping `tp:A:P` (primary/supplementary) and `tp:A:I`
  (inversion), dropping `tp:A:S`. One hit per line.
- **Alignment file** (`src/aln.rs::run`): the container (SAM/BAM/CRAM) is
  autodetected by `noodles_util::alignment::io::reader` (magic-byte sniff, works
  on stdin too); records come as `Box<dyn sam::alignment::Record>` and are
  processed by the trait (one path for all three). Iterate primary records only
  (skip secondary/supplementary/unmapped). Each read's hits = the primary (from
  its own CIGAR) + one per `SA:Z:` entry. For paired-end reads the qname gets a
  `/1` (first segment) or `/2` (last segment) suffix so the two mates are
  distinct; unpaired reads keep the bare name. No grouping needed, so sorted and
  name-grouped files both work. Coordinates are computed to match PAF exactly
  (`span_from_ops`): `qspan/refspan/alen` from the CIGAR, `lead/tail` clips,
  `qlen=lead+qspan+tail`, `qs = +?lead:tail`, `qe=qs+qspan`, `ts=pos0`,
  `te=pos0+refspan`.
  - **CRAM** needs a reference: `-r` (a faidx-indexed FASTA, i.e. `ref.fa` +
    `ref.fa.fai`) builds a `fasta::Repository` the reader decodes sequences
    against. A 4-byte `CRAM` magic peek gates this: CRAM without `-r` errors; `-r`
    is ignored (repository not built) for BAM/SAM. `eseq` for CRAM is the
    reference-reconstructed `SEQ`, identical to the BAM output for the same reads.
  - **ALT contigs** (`--alt FILE`, alignment-file input only; `load_alt`): a
    `HashSet` of contig names — `@`-prefixed lines (SAM headers) are ignored and
    the first tab-delimited column of the rest is the contig name, so a plain
    one-name-per-line list or a bwa-kit `.alt` SAM both work (gzip ok). A read
    whose **primary** hit lands on an ALT contig is skipped entirely; ALT hits in
    the `SA` tag are dropped before clips/joins form. No-op for `--paf`.
  - noodles crates: util 0.82 (feature `alignment`) / sam 0.87 / core 0.20 /
    cram 0.96 / fasta 0.64.

Within a read, hits are sorted by query start `qs`. Filtering is off by default.
`-a` drops hits whose alignment block length (PAF col 11, the `alen` field) is
below a threshold (default 0) — the one **upfront** hit filter (`passes_filter`),
applied before breakends are formed. The clip threshold is `-c` (default 50).
`-s` sets the `source=` dataset label (default `foo`), stamped on every record
(alignment file and PAF).

Two **post-filters** drop emitted lines (never alignments) — like their `merge`
namesakes — so the per-read `idx`/`n_aln` counters are unchanged (dropped
breakends are simply absent, so `idx` stays a stable per-breakend identifier):
`-Q` (default 0 = keep all) drops breakends whose computed `equal` quality is
below the threshold (breakends without an `equal` value — no `eseq`, no `QUAL`,
or PAF — are always kept); `-q` (default 0) drops lines whose col-7 mapq (the
join's `max`, or a clip's mapq) is below the threshold. Both shrink the output
feeding downstream tools. In the `--help` listing `-q` is placed after `-Q` to
reflect that it is a post-filter.

**End-of-run stats** (`extract.rs::Stats`, printed to stderr after a successful
run): number of reads, total bases in primary alignments, and read N50. A read
is one primary record (alignment input; counted after the ALT skip) or one PAF
read group; counted **before** the `-a` filter. Primary bases sum the primary
hit's aligned query span (`qe - qs`); PAF carries no primary flag (all kept
lines are `tp:A:P`), so the read's longest query span stands in. N50: walking
reads from longest to shortest, the read length at which the running sum first
reaches half the total read length. To stay O(1)-memory on billion-read
short-read BAMs, lengths `<= 65535` (`MAX_COUNTED_LEN`) are tallied in a fixed
counting array and only longer reads are kept individually; the N50 walk visits
the long list first, then the array from the top (unit test:
`extract::tests::n50_matches_naive`).

**eseq/elen (alignment-file input only, not PAF).** With the read sequence
available, each breakend also gets a window around the junction: read-forward
`[max(lo-f,0), min(hi+f,qlen)]` where `lo=min(y0.qe,y1.qs)`, `hi=max(..)` for a
join (`lo=hi=p`, the clip point `first.qs`/`last.qe`, for a clip). `-f` (default
250) is the flank, `-e` (default 1000) the max window.
`elen=leftFlank,qdist,rightFlank` is always emitted; `eseq=<bases>` is appended
only when the window `<= -e`. Both are read-forward (NOT flipped with the output,
unlike `qlen`). `eseq` is on the original read strand: the primary `SEQ`
reverse-complemented iff the primary is reverse-mapped (`src/aln.rs::revcomp`),
usable only when `SEQ.len()==qlen` (soft-clipped primary). For CRAM the `SEQ` is
reconstructed against the `-r` reference. See `src/extract.rs::eseq_info`.

`equal=<phred>` is the quality of the eseq window, emitted **between `elen` and
`eseq`** — only when `eseq` is emitted and per-base qualities cover the window
(omitted when `QUAL` is `*`). It turns each base quality to an error rate, means
the expected errors over the window, and scales back to phred:
`round(-10*log10((Σ 10^(-Q_i/10))/L))`, `L=|eseq|` (`eseq_qual`). Qualities are
read-forward too: reversed (not complemented) on a reverse-mapped primary. For
constant quality `Q` this is `Q`.

### Alignment file vs PAF parity

For the same alignments the output is identical except `aln_len` of supplementary
(SA-derived) hits, which is a few bp smaller from a BAM/CRAM because minimap2's
`SA` CIGAR merges small indels (the primary hit's `aln_len` matches PAF).
Everything else — line count, offsets, `ori`, `strand`, `qlen`, `mapq`, clip/join
logic, and inversion (`tp:A:I`) segments — matches exactly. (Verified on the test
data: both emit 1009 lines; the diff after masking `aln_len` is empty.) BAM and
CRAM outputs are byte-identical for the same reads (`tests/extract.rs::cram_matches_bam`).

### Output (9 columns, TAB-delimited)

```
ctg1  pos1  ori  ctg2  pos2  qname  mapq  strand  source=..;idx=..;n_aln=..;aln_len=..;qlen=..;mapq=..[;elen=..[;equal=..;eseq=..]]
```

`ori` is two characters, each `>`/`<`, or `.` for a missing mate. Tag **order in
the INFO column is not significant** — every consumer parses by key, not
position — but the emission order is `source`, `idx`, then three tags keyed to
**output order** — index 1 tracks `ctg1:pos1`, index 2 tracks `ctg2:pos2` — so
their values swap with the endpoints when a join is flipped (`source` is
per-read and `elen`/`eseq`, appended for BAM, do not swap):

- `source=NAME` — the `-s` dataset label (default `foo`), for telling reads of
  different datasets apart (tumor/normal, trio) downstream; written first purely
  for readability. Same on every record of a run.
- `idx=N` — 0-based index of this clip/breakend within the read, emission order
  (left clip, joins, right clip); resets per read (`emit_read`'s local counter).
- `n_aln=N` — number of alignments (primary + supplementary) for the read that
  passed the `-q`/`-a` filters, i.e. `hits.len()` at `emit_read`. A read-level
  constant, the same on every breakend of the read.
- `aln_len=len1,len2` — alignment block length (PAF col 11) of each hit; for a
  clip `len2 = 0`.
- `qlen=left,middle,right` — query lengths summing to the read length.
  `left` is on the `ctg1:pos1` side, `right` on the `ctg2:pos2` side. Clip:
  `right` = clipped length (`qs` for a left clip, `qlen-qe` for a right clip),
  `middle` = 0. Join: `left = y0.qe`, `middle = y1.qs - y0.qe` (query gap, may be
  negative), `right = qlen - y1.qs`, computed in read order then swapped on flip.
- `mapq=mapq1,mapq2` — mapq of each hit; for a clip `mapq2 = 0`. Note this is
  distinct from the field-7 mapq column, which is `max(y0.mapq, y1.mapq)` for a
  join.

- **Left clip** (`first.qs > c`): `first.ctg fts(first) arrow_l.  .  .  qname mapq strand INFO`
- **Right clip** (`last.qlen - last.qe > c`): `last.ctg fte(last) arrow_r.  .  .  qname mapq strand INFO`
- **Join** (each adjacent pair `y0`→`y1`): endpoint `c0 = (y0.ctg, fte(y0), ori(y0), y0.alen)`
  joins `c1 = (y1.ctg, fts(y1), ori(y1), y1.alen)`. Canonicalize so the smaller
  `(ctg, pos)` is first; if swapped, invert both arrows, swap the two `alen`
  values, the two `mapq` values, and the `qlen` `left`/`right`, and set the
  `strand` column to `-` (else `+`). The field-7 `mapq` column is `max(y0, y1)`.

Helpers: `fts = strand=='+'? ts : te`, `fte = strand=='+'? te : ts`,
`ori = strand=='+'? '>' : '<'`. Clip arrows: `arrow_l = flip(ori)`,
`arrow_r = ori`.

Emission order per read: left clip, joins (read order), right clip.

### Note on `.msv` fixtures

The `test/join*.msv` files carry an INFO column
(`SVTYPE=BND;...;aln_len=...;source=foo`) produced by `minisv.js` that differs
from badclip's — e.g. minisv's `aln_len` uses query-span in read order, badclip
uses PAF col-11 block length in output order. So tests compare against each
`.msv` truncated to its first 8 columns with badclip's own INFO column
(`aln_len=...;qlen=...;mapq=...`) appended (see `tests/extract.rs`). BAM fixtures
(`test/bam01.bam`, `test/bam01.srt.bam`, three reads) are compared against the
committed golden output `test/bam01.expected` (eseq bytes verified against
samtools); `test/bam01.paf` is the single-alignment read used for BAM/PAF parity;
`test/inv01.paf` is a `tp:A:I` inversion read (two joins).

## Not implemented (deliberately, for now)

- CIGAR-based indel (INS/DEL) extraction.
- Same-contig SV typing (INS/DEL/DUP/INV) and the trailing INFO column — always
  emits `BND`-style orientation regardless of contig.
