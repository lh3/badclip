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
- `src/bam.rs`    — BAM input via noodles: build `Hit`s from a primary record + its `SA` tag.
- `src/extract.rs`— dispatch (BAM vs `--paf`), grouping, sorting, clip/join emission; shared `emit_read`/`passes_filter`.
- `tests/extract.rs` — end-to-end tests driving the compiled binary.
- `test/`         — `*.paf` inputs, `*.msv` (minisv) references, `bam01*.bam` fixtures, `minisv.js`.

## Interval / offset notation

An interval is a pair of **offsets** `|st,en|`, where an offset sits *between*
two adjacent bases, not on a base (0-based, half-open). In `"A|CGT|AGC"`,
`|1,4|` is `"CGT"`. A hit matches a target interval `tc:|ts,te|` to a query
interval `|qs,qe|` on a strand. Emitted positions are raw 0-based offsets — no
±1 adjustment.

## `extract` — behaviour

Input is **BAM by default**, or **PAF with `--paf`**. `-` reads stdin (no input
at all → print help, don't block on stdin). Both paths build a `Vec<Hit>` per
read and feed the shared `emit_read`.

- **PAF** (`src/extract.rs::run_paf`): assumed grouped by read name; streamed and
  buffered per read, keeping `tp:A:P` (primary/supplementary) and `tp:A:I`
  (inversion), dropping `tp:A:S`. One hit per line.
- **BAM** (`src/bam.rs::run_bam`): iterate primary records only (skip
  secondary/supplementary/unmapped). Each read's hits = the primary (from its
  own CIGAR) + one per `SA:Z:` entry. No grouping needed, so sorted and
  name-grouped BAM both work. Coordinates are computed to match PAF exactly
  (`span_from_ops`): `qspan/refspan/alen` from the CIGAR, `lead/tail` clips,
  `qlen=lead+qspan+tail`, `qs = +?lead:tail`, `qe=qs+qspan`, `ts=pos0`,
  `te=pos0+refspan`. noodles crates: bam 0.92 / sam 0.87 / core 0.20.

Within a read, hits are sorted by query start `qs`. Filtering is off by default:
`-q` drops hits below a mapq (default 0) and `-a` drops hits whose alignment
block length (PAF col 11, the `alen` field) is below a threshold (default 0).
The clip threshold is `-c` (default 100).

**eseq/elen (BAM only).** With the read sequence available, each breakend also
gets a window around the junction: read-forward `[max(lo-f,0), min(hi+f,qlen)]`
where `lo=min(y0.qe,y1.qs)`, `hi=max(..)` for a join (`lo=hi=p`, the clip point
`first.qs`/`last.qe`, for a clip). `-f` (default 250) is the flank, `-e` (default
1000) the max window. `elen=leftFlank,qdist,rightFlank` is always emitted (BAM);
`eseq=<bases>` is appended only when the window `<= -e`. Both are read-forward
(NOT flipped with the output, unlike `qlen`). `eseq` is on the original read
strand: the primary `SEQ` reverse-complemented iff the primary is reverse-mapped
(`src/bam.rs::revcomp`), usable only when `SEQ.len()==qlen` (soft-clipped
primary). See `src/extract.rs::eseq_info`.

### BAM vs PAF parity

For the same alignments the output is identical except `aln_len` of supplementary
(SA-derived) hits, which is a few bp smaller from BAM because minimap2's `SA`
CIGAR merges small indels (the primary hit's `aln_len` matches PAF). Everything
else — line count, offsets, `ori`, `strand`, `qlen`, `mapq`, clip/join logic,
and inversion (`tp:A:I`) segments — matches exactly. (Verified on the test data:
both emit 1009 lines; the diff after masking `aln_len` is empty.)

### Output (9 columns, TAB-delimited)

```
ctg1  pos1  ori  ctg2  pos2  qname  mapq  strand  aln_len=..;qlen=..;mapq=..[;elen=..;eseq=..]
```

`ori` is two characters, each `>`/`<`, or `.` for a missing mate. The INFO
column's first three tags are keyed to **output order** — index 1 tracks
`ctg1:pos1`, index 2 tracks `ctg2:pos2` — so their values swap with the endpoints
when a join is flipped (`elen`/`eseq`, appended for BAM, do not swap):

- `aln_len=len1,len2` — alignment block length (PAF col 11) of each hit; for a
  clip `len2 = 0`.
- `qlen=left,middle,right` — query lengths summing to the read length.
  `left` is on the `ctg1:pos1` side, `right` on the `ctg2:pos2` side. Clip:
  `right` = clipped length (`qs` for a left clip, `qlen-qe` for a right clip),
  `middle` = 0. Join: `left = y0.qe`, `middle = y1.qs - y0.qe` (query gap, may be
  negative), `right = qlen - y1.qs`, computed in read order then swapped on flip.
- `mapq=mapq1,mapq2` — mapq of each hit; for a clip `mapq2 = 0`. Note this is
  distinct from the field-7 mapq column, which is `min(y0.mapq, y1.mapq)` for a
  join.

- **Left clip** (`first.qs > c`): `first.ctg fts(first) arrow_l.  .  .  qname mapq strand INFO`
- **Right clip** (`last.qlen - last.qe > c`): `last.ctg fte(last) arrow_r.  .  .  qname mapq strand INFO`
- **Join** (each adjacent pair `y0`→`y1`): endpoint `c0 = (y0.ctg, fte(y0), ori(y0), y0.alen)`
  joins `c1 = (y1.ctg, fts(y1), ori(y1), y1.alen)`. Canonicalize so the smaller
  `(ctg, pos)` is first; if swapped, invert both arrows, swap the two `alen`
  values, the two `mapq` values, and the `qlen` `left`/`right`, and set the
  `strand` column to `-` (else `+`). The field-7 `mapq` column is `min(y0, y1)`.

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
