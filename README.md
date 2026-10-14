> [!Warning]
> This project is vibe coded with Claude Code.

# badclip

Extract breakends and structural-variant signals from long-read alignments.

`badclip` reads read-to-reference alignments and reports **breakends**: the
junctions between chimeric alignments of the same read (*splits*), long
deletions and insertions inside an alignment's CIGAR (*del*/*ins*), and,
optionally, read ends that are soft-clipped away from every alignment
(*trims*). Splits reveal structural rearrangements (translocations, inversions,
large indels); trims flag reads that end abruptly against nothing.

It has six subcommands, forming a pipeline:

| Subcommand | Purpose |
|------------|---------|
| [`extract`](#extract) | Find breakends in read alignments (the main step). |
| [`getfa`](#getfa) | Turn `extract` output into a FASTA of breakend sequences. |
| [`fltbnd`](#fltbnd) | Drop breakends already explained by a pangenome. |
| [`merge`](#merge) | Cluster per-read breakends into consensus SV calls. |
| [`fltreg`](#fltreg) | Drop breakends that fall in BED regions. |
| [`fltcnt`](#fltcnt) | Keep `merge` calls with enough raw/filtered read support. |

Every subcommand reads its input from a file or from stdin via `-`, transparently
decompresses gzip'd input, and prints its own `--help` (instead of waiting on
stdin) when run with no input. `getfa` and `fltbnd` were called `geteseq` and
`flteseq` before v0.2.0; the old names still work as hidden aliases.

## Install

```sh
cargo build --release
# binary at target/release/badclip
```

## `extract`

Find breakends (splits, CIGAR deletions/insertions and, with `-T`, trims) in
read-to-reference alignments.

```sh
badclip extract [OPTIONS] [INPUT]
```

- `INPUT` — an **alignment file** (SAM/BAM/CRAM, format auto-detected) by
  default, or **PAF** with `--paf`. Pass `-` to read from stdin (autodetection
  works on a pipe too). PAF input may be gzip'd (auto-detected). **CRAM requires a
  reference** via `-r`.

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `--paf` | off | Read PAF (optionally gzip'd) instead of an alignment file. |
| `-r`, `--reference <FASTA>` | — | Reference FASTA (faidx-indexed), **required for CRAM**; ignored for BAM/SAM/PAF. |
| `--alt <FILE>` | — | File of ALT contig names (alignment-file input): `@` header lines are ignored, the first tab-delimited column of the rest is the contig name (a plain list or a bwa-kit `.alt` SAM both work). Reads whose primary maps to an ALT contig are skipped, and ALT hits are dropped before breakends are formed. |
| `-s`, `--source <STR>` | `foo` | Dataset name, stamped on every record as a `source=` INFO tag. |
| `-T`, `--trim` | off | Also report **trim** breakends (read ends clipped away from every alignment). |
| `-c`, `--min-clip <INT>` | `50` | Minimum clip length to report a trim breakend. |
| `-g`, `--min-gap <NUM>` | `1000` | Minimum CIGAR deletion/insertion length to report a del/ins breakend (`k`/`m`/`g` suffix ok); `0` disables. |
| `-a`, `--min-aln-len <INT>` | `0` | Drop hits whose alignment block length (PAF col 11) is below this value (upfront hit filter). |
| `-f`, `--flank <INT>` | `250` | Flanking read sequence extracted on each side of a breakend (alignment-file input only). |
| `-e`, `--max-eseq <INT>` | `1000` | Maximum extracted window; longer windows omit `eseq` (alignment-file input only). |
| `-Q`, `--min-equal <INT>` | `0` | **Post-filter:** drop breakends whose eseq quality (`equal`; for an insertion the lower of `equal`/`equal2`) is below this, to shrink downstream input. `0` keeps all; breakends without an `equal` value are always kept. |
| `-q`, `--min-mapq <INT>` | `0` | **Post-filter:** drop output lines whose col-7 mapq (the split's larger, or the one hit's mapq for trim/del/ins) is below this. Does not drop alignments upfront. |

Examples:

```sh
badclip extract aln.bam                     # BAM (default), sorted or unsorted
badclip extract -T aln.bam                  # also report trims
badclip extract aln.cram -r ref.fa          # CRAM (needs faidx-indexed reference)
samtools view -b … | badclip extract -      # BAM/CRAM from stdin (autodetected)
badclip extract --paf aln.paf.gz            # PAF (gzip auto-detected)
minimap2 -c … | badclip extract --paf -     # PAF from stdin (-c: cg:Z for del/ins)
```

At the end of a run, `extract` prints input statistics to stderr: the number of
reads (primary records, including unmapped ones; secondary and supplementary
records are ignored), the number of unmapped reads (for PAF, lines whose target
name is `*`, as from `minimap2 --paf-no-hit`), the total bases in primary
alignments (the primary hit's aligned query span per read; for PAF, which has no
primary flag, the read's longest span), and the read N50 length. Reads are counted before the `-a` filter.

### Input assumptions

- **Alignment file (SAM/BAM/CRAM)**: the container format is auto-detected. A
  read's chimeric hits (for splits and trims) are read from the primary
  alignment's `SA:Z:` tag. The `SA` CIGARs are collapsed, so del/ins breakends
  come from each record's own CIGAR — the primary's and every supplementary
  record's (a supplementary record contributes only its del/ins). Secondary and
  unmapped records are ignored. Each record is handled on its own, so this works
  on both coordinate-sorted and name-grouped files (no grouping required).
  Paired-end reads get a `/1` or `/2` suffix on the read name (first/last
  segment) so the mates stay distinct; unpaired reads keep the bare name.
  **CRAM** needs `-r ref.fa` (with `ref.fa.fai`) to decode read sequences; `-r`
  is ignored for BAM/SAM.
- **PAF**: alignments are grouped by read name (as produced directly by
  minimap2); primary/supplementary (`tp:A:P`) and inversion (`tp:A:I`)
  alignments are used, secondary (`tp:A:S`) are ignored. Del/ins come from the
  `cg:Z:` tag (minimap2 `-c`); without it there are none.
- Within a read, hits are sorted by their start position along the read.

An alignment file and PAF produce the same breakends for the same alignments,
with two caveats: `aln_len` can differ (supplementary hits are a few bp smaller
from BAM/CRAM because minimap2's `SA` CIGAR is collapsed, and minimap2's PAF
column 11 can differ from the CIGAR on alignments with long indels), and PAF
insertions are not normalized (no read sequence; see below), so they sit at the
CIGAR position.

### Output

Tab-delimited, 9 columns:

```
ctg1   pos1   ori   ctg2   pos2   qname   mapq   strand   source=..;ctype=..;idx=..;n_aln=..;aln_len=..;qlen=..;mapq=..[;inslen=..][;elen=..[;equal=..;eseq=..]][;elen2=..[;equal2=..;eseq2=..]]
```

`ori` is two characters, each `>` or `<` (or `.` for a missing mate).
Positions are raw 0-based offsets that sit *between* bases. The INFO column is a
list of `key=value` tags; **their order is not significant** (consumers parse by
key). It is written as `source`, `ctype`, `idx`, `n_aln`, then three tags in
output order (index 1 ↔ `ctg1:pos1`, index 2 ↔ `ctg2:pos2`):

- `source=NAME` — the `-s` dataset label (default `foo`), for telling reads from
  different datasets apart (tumor/normal, trio) in downstream tools. Written
  first only for readability; the same on every record.
- `ctype=TYPE` — the breakend type: `trim` (a read end not aligned; only with
  `-T`), `split` (a junction between two alignments), `del` or `ins` (a CIGAR
  deletion/insertion of at least `-g` bases).
- `idx=N` — identifies the breakend within the read. For trims/splits it is the
  0-based index in emission order (left trim, splits, right trim), resetting per
  read; a hidden trim (no `-T`) still takes its number. For del/ins it is the
  event's read-forward query offset (a read's del/ins can come from several
  records, so they cannot share that counter).
- `n_aln=N` — number of alignments (primary + supplementary) for the read that
  passed the `-a` filter. Same on every breakend of the read (e.g. `1` for a
  single-alignment read, `3` for a 3-way chimera).
- `aln_len=len1,len2` — alignment length (PAF column 11) of each hit; a trim has
  no second hit, so `len2 = 0`; for del/ins both are the one alignment's.
- `qlen=left,middle,right` — query lengths that sum to the read length. For a
  trim, `right` is the clipped length and `middle` is `0`; for a split, `middle`
  is the (possibly negative) query gap between the two hits; for an insertion,
  `middle` is the inserted length (`0` for a deletion).
- `mapq=mapq1,mapq2` — mapq of each hit; for a trim `mapq2 = 0`. (The `mapq`
  *column* above is the larger of a split's two mapqs; this tag keeps both.)
- `inslen=N` — (insertions only) the inserted length.

The remaining tags come from **alignment-file input only** (PAF has no read
sequence; for CRAM the sequence is reconstructed against `-r`), in read-forward
orientation (not flipped with the output):

- `elen=leftFlank,qdist,rightFlank` — sizes of the extracted window: up to `-f`
  bases either side of the junction (`qdist` is the query gap of a split, `0`
  otherwise), so `|eseq| = leftFlank + |qdist| + rightFlank`. The window is
  clamped to the bases the record carries (a hard-clipped supplementary record
  only has its aligned part).
- `equal=<phred>` — the quality of the `eseq` window: each base quality is turned
  into an error rate, the expected errors are averaged over the window and scaled
  back to phred (`-10·log10((Σ 10^(-Q/10))/|eseq|)`). Sits between `elen` and
  `eseq`, emitted only when `eseq` is (and when the read has base qualities). For
  a constant quality `Q` it is `Q`.
- `eseq=<bases>` — the read sequence over that window, on the **original read
  strand**. Omitted (while `elen` stays) when the window exceeds `-e`.
- `elen2`, `equal2`, `eseq2` — (insertions only) the same for the insertion's
  second breakpoint on the read; `elen`/`equal`/`eseq` cover the first
  (read-forward order: on a reverse-strand alignment the first is the
  reference-right junction).

- **Trim** (`-T`) — one side is a real breakend, the other is nothing (`.`):

  ```
  chr1   123569841   >.   .   .   read/ccs   1   -   source=foo;ctype=trim;idx=0;aln_len=14732,0;qlen=14925,0,141;mapq=1,0
  chr1   123555134   <.   .   .   read/ccs   1   -   source=foo;ctype=trim;idx=1;aln_len=14732,0;qlen=14774,0,292;mapq=1,0
  ```

- **Split** — two mapped ends meet. The smaller `(contig, position)` is written
  first; the two arrows show how the sides are oriented:

  ```
  chr1    57375269   >>   chr21   32069271   read/ccs   60   -   source=foo;ctype=split;idx=0;aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60
  chr13   51911798   <>   chr2    76026062   read/ccs   60   +   source=foo;ctype=split;idx=0;aln_len=29505,5668;qlen=29436,1,5661;mapq=60,60
  ```

  `chr1:57375269 >> chr21:32069271` means the **right** side of `chr1:57375269`
  joins the **left** side of `chr21:32069271`. `chr13:51911798 <> chr2:76026062`
  means the **left** side of `chr13:51911798` joins the **left** side of
  `chr2:76026062`.

- **Del** — a CIGAR deletion of `|t,t+len|`, written like the split an aligner
  would report for it (`>>`, column 8 = the alignment strand):

  ```
  chr7   31550975   >>   chr7   31551986   read/ccs   60   +   source=foo;ctype=del;idx=23562;n_aln=1;aln_len=34703,34703;qlen=23562,0,9633;mapq=60,60;elen=250,0,250;equal=25;eseq=…
  ```

- **Ins** — a CIGAR insertion, written `ctg L >> ctg R`, where `L`/`R` are the
  leftmost/rightmost reference positions it can be placed at equivalently: the
  insertion is shifted right while its first base equals the read base after it
  (then its second base, …), and left likewise. For an insertion with a target
  site duplication (e.g. a mobile element), `R - L` is about the TSD length.
  Shifts stay within the matched bases adjacent to the insertion (stopping at
  another indel or an `N`) and use the read sequence only, so PAF insertions
  are not shifted (`L = R`). Deletions are not shifted either (that would need
  the deleted reference bases).

  ```
  chr7   93405377   >>   chr7   93405381   read/ccs   60   +   source=foo;ctype=ins;idx=8478;n_aln=1;aln_len=19742,19742;qlen=8478,1802,9449;mapq=60,60;inslen=1802;elen=250,0,250;equal=28;eseq=…;elen2=250,0,250;equal2=28;eseq2=…
  ```

## `getfa`

Convert `extract` output into a FASTA of the extracted breakend sequences
(formerly `geteseq`, still accepted as an alias).

```sh
badclip getfa [OPTIONS] [INPUT]
```

- `INPUT` — `extract` output (gzip ok; `-` or omit for stdin).

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `-Q`, `--min-equal <INT>` | `20` | Drop records whose eseq quality (`equal` tag; for an insertion the lower of `equal`/`equal2`) is below this. Records without an `equal` tag are kept. |

For each eseq window of a record, it writes one FASTA entry named
`readName_idx_leftFlank_rightFlank_code` (from the `idx` and `elen` tags) with
the window's bases; records without an eseq are skipped. `code` is the breakend
type: `T` (trim), `S` (split), `D` (del), or, for an insertion's two windows,
`L` (`eseq`, the first breakpoint along the read) and `R` (`eseq2`, the
second). Records without a `ctype` tag (pre-0.2 output) count as a trim when
they have no mate (column 4 = `.`), else as a split.

```sh
badclip extract aln.bam | badclip getfa - > eseq.fa   # "-" = stdin
badclip getfa extract.txt                             # or a file
```

Example output:

```
>m84039_.../234884533/ccs_0_250_250_T
ACTTTGGGAGGCCAAGGCAGGCGGATCACCTG...
```

## `fltbnd`

Keep only breakends whose junction sequence is **not** already explained by a
pangenome (formerly `flteseq`, still accepted as an alias). Align the `getfa`
FASTA against a pangenome index (e.g. `ropebwt3 sw`) to get a PAF, then:

```sh
badclip fltbnd [OPTIONS] <extractOut> <ropebwt3.paf>
```

- `<extractOut>` — `extract` output (gzip ok; `-` for stdin).
- `<ropebwt3.paf>` — the `ropebwt3 sw` PAF of the `getfa` FASTA (gzip ok).

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `-l`, `--margin <INT>` | `50` | Margin extending the protected junction interval on each side. |
| `-s`, `--source <STR>` | off | Print **all** input lines and rewrite the `source=` tag to `STR` on survivors (kept/novel breakends). |
| `-Q`, `--min-equal <INT>` | `20` | Drop lines whose eseq quality (`equal` tag; for an insertion the lower of `equal`/`equal2`) is below this (both modes). Lines without an `equal` tag are kept. |

A record is **dropped** if it has no `eseq` tag, or if some ropebwt3 alignment's
query interval `[qs,qe]` on an eseq window contains that window's junction
interval `[max(0, elen0-l), min(elen0+|elen1|+l, |eseq|)]` (i.e. the alignment
spans the breakend, with `-l` margin on each side — the breakend is
"protected"). An insertion has two windows (`_L`, `_R`) and is dropped if
**either** junction is spanned. Surviving `extract` lines are printed verbatim.

With `-s STR`, the output instead contains **every** input line, with the
survivors' `source=` tag rewritten to `STR` (dropped lines keep their original
`source`). This marks novel breakends distinctly so a downstream `merge` can tell
them from pangenome-explained ones.

The PAF's query names must be the `getfa` names, and in the same order as the
FASTA records (as produced by the pipeline), so the two files are streamed
without loading into memory.

```sh
badclip extract aln.bam > clip.txt
badclip getfa clip.txt | ropebwt3 sw -d pangenome.fmd - > rb3.paf
badclip fltbnd clip.txt rb3.paf > novel.txt
```

## `merge`

Cluster per-read breakends from `extract` into consensus SV calls, each with a
supporting-read count. Unlike the streaming subcommands above, `merge` loads its
input into memory and sorts it by `(ctg1, pos1)`, so the input need not be
pre-sorted.

```sh
badclip merge [OPTIONS] [INPUT]...
```

- `INPUT` — `extract` output (gzip ok; `-` for stdin). Several files are merged
  as if concatenated, e.g. `badclip merge a.clip.gz <(badclip fltbnd ...)`.

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `-c`, `--min-cnt <INT>` | `2` | Minimum supporting-read count to emit a call. |
| `-s`, `--min-cnt-strand <INT>` | `0` | Minimum read count on each strand. |
| `-w`, `--win-size <INT>` | `100` | Clustering window (bp). |
| `-A`, `--max-allele <INT>` | `100` | Cap on simultaneously-open clusters. |
| `-M`, `--max-check <INT>` | `500` | Maximum reads compared per cluster. |
| `-Q`, `--min-equal <INT>` | `20` | Drop input breakends whose eseq quality (`equal` tag; for an insertion the lower of `equal`/`equal2`) is below this before clustering. Breakends without an `equal` tag are kept. |
| `-q`, `--min-mapq <INT>` | `0` | Drop input breakends whose col-7 mapq (the split's larger, or a trim's mapq) is below this. |
| `-p`, `--min-mapq-min <INT>` | `10` | Drop input breakends with a mate (splits, del/ins) whose *smaller* per-hit mapq (the `min` of the `mapq=` pair) is below this; trims are exempt. |
| `-m`, `--merge-input` | off | Sample-merge mode: treat input as combined `merge` output (see below). |
| `-C`, `--min-cnt-in <INT>` | `3` | (`-m` only) drop input lines whose total count is below this. |
| `-S`, `--min-cnt-strand-in <INT>` | `1` | (`-m` only) drop input lines whose per-strand count is below this. |

Breakends within `-w` bp on both endpoints (with `><`/`<>` treated as the same
inversion) are grouped, whatever their `ctype` — a CIGAR deletion clusters with
the split-read deletions at the same coordinates, and the call reports the mix in
`ctype_count`; a group is reported only if it has at least `-c` reads
with at least `-s` on each strand. Output keeps the 9-column, TAB-delimited
`extract` shape, with the supporting-read count in the `mapq` slot and a
merge-derived INFO column:

```
ctg1   pos1   ori   ctg2   pos2   .   count   strand   avg_mapq=q1,q2;count=..;ctype_count=..[;avg_inslen=..;diff_inslen=..][;count_fr=..;count_rf=..[;foldback]][;ins_nearby=..];reads=..
```

- `avg_mapq=q1,q2` — mean per-hit mapq of the cluster's reads, split by output
  endpoint: `q1` on the `ctg1:pos1` side, `q2` on the `ctg2:pos2` side (`q2 = 0`
  for a trim cluster).
- `count=src:fwd,rev|…` — for each `source=` dataset present in the cluster (from
  the reads' `source=` tag), its forward,reverse read counts; `|`-joined and
  sorted alphabetically by source. Only sources actually present are listed
  (never `x:0,0`), e.g. `count=foo:2,3|retain:4,5`, or `count=foo:2,2` for a
  single source. The column-7 `count` is the grand total across all sources and
  strands.
- `ctype_count=type:n|…` — the reads per breakend type (`del`, `ins`, `split`,
  `trim`), alphabetical, listing only types present, e.g. `del:10|split:5` for a
  deletion seen both in CIGARs and as split reads. Sums to column 7. Input
  without a `ctype` tag (pre-0.2) counts as `trim` when it has no mate, else
  `split`.
- `avg_inslen`, `diff_inslen` — for a call with `ins` reads: their mean
  `inslen` (rounded) and `max - min`.
- `count_fr`, `count_rf` — reads oriented `><` and `<>`; present only when the
  cluster has at least one such read. `foldback` is added when only one of the
  two orientations is present (`count_fr` or `count_rf` is `0`) and both
  endpoints are on the same contig.
- `ins_nearby=src:n|…` — for a call with `split`/`trim` reads: the number of
  distinct `ins` reads in the whole input (any source, after `-Q`/`-q`/`-p`,
  whether or not their own call is reported) within `-w` bp of either printed
  end — an `ins` line `ctg L >> ctg R` is near `(ctg, p)` when
  `L - w <= p <= R + w` — per source, alphabetical; omitted when none. The
  call's own reads (same source and name) are not counted, so an insertion
  seen as both `split` and `ins` reads does not flag itself. A long
  germline insertion whose sequence maps elsewhere can masquerade as a
  translocation; `ins` reads nearby (e.g. from the matched normal) flag it.
- `reads=src:name,…|…` — the supporting read names, stratified by `source=` the
  same way as `count=`: `|`-joined and alphabetical by source (names in
  member order within each source), e.g. `reads=foo:r1,r2|retain:r3,r4`. Lists the
  same sources, in the same order, as `count=`.

```sh
badclip extract aln.bam | badclip merge - > sv.txt
```

### `-m` sample-merge

With `-m`, `merge` re-runs the same clustering on the (filtered) combination of
several samples' `merge` outputs — a second-level, population merge. Only the
INFO parsing and emission change; the clustering is identical. mapq is read from
`avg_mapq=q1,q2` and per-strand read counts from the `count=` tag (summed over
sources); the `-c`/`-s` cluster filters and `-q`/`-p` stay read-based (`-Q` is a
no-op — merge output has no `equal` tag). The per-line `-C`/`-S` filters drop
under-supported input lines before clustering.

Output stays in the 9-column shape but with **column 7 = the number of distinct
sources** (samples) in the cluster and a compact INFO:
`avg_mapq=q1,q2;count=F,R` (aggregate forward,reverse reads, no per-source
breakdown) plus `count_fr`/`count_rf`/`foldback` for inversions — **no `reads=`
tag**, and none of `ctype_count`/`avg_inslen`/`diff_inslen`/`ins_nearby`
(ignored on input). The output is itself `-m`-parseable, so merges chain.

```sh
# combine per-sample merge outputs, then re-merge across samples
cat sample*.merge.txt | badclip merge -m - > population.txt
```

## `fltreg`

Drop `extract`/`merge` breakends that fall in BED regions (e.g. blacklists,
centromeres, low-complexity regions).

```sh
badclip fltreg [-l INT] [INPUT] [BED]
```

- `INPUT` — `extract` or `merge` output (gzip ok; `-` or omit for stdin).
- `BED` — BED file of regions to filter against (gzip ok).
- `-l INT` — pad each BED interval to `[start-l, end+l)` before testing
  (default 0, i.e. the raw intervals).

Both output formats share the 9-column layout, so one filter serves both. A line
is **dropped** if **either** breakend endpoint — `(ctg1, pos1)` (columns 1,2) or,
when present, `(ctg2, pos2)` (columns 4,5) — falls in a region; a trim
(`ctg2 = "."`) tests only its one endpoint. Survivors are printed verbatim; a
line whose endpoints can't be parsed is kept.

Positions are raw 0-based offsets and BED intervals are half-open, so an offset
`p` is inside `[start, end)` iff `start <= p < end` (or `start-l <= p < end+l`
with `-l`). The BED is loaded into per-contig intervals, padded by `-l`, sorted
by start and merged into a disjoint list, so each
lookup is a single binary search (no new dependency; overlapping BED input is
handled).

```sh
badclip merge - < sv.txt | badclip fltreg - blacklist.bed.gz > sv.filtered.txt
```

## `fltcnt`

Keep `merge` calls with enough read support from selected sources — by default
the `fltbnd`-filtered (`*.flt`) ones.

```sh
badclip fltcnt [OPTIONS] [INPUT]
```

- `INPUT` — `merge` output (gzip ok; `-` or omit for stdin).
- `-c INT` — minimum reads from the `--src` sources (default 3).
- `-s INT` — minimum `--src` reads on each strand (default 0).
- `-r INT` — maximum reads from the `--rest` sources (default 0, i.e. every
  read must come from `--src`; `-1` = no limit).
- `-l NUM` — minimum `pos2 - pos1` (columns 5 and 2). Any value `>= 0` also
  drops trims (column 4 = `.`), so the default 0 drops trims only; `-1` disables
  the filter. Inter-contig splits have no distance and pass. Note an insertion
  spans only its `L..R` shift range, so any positive `-l` drops insertions.
  Accepts a `k`/`m`/`g` suffix, e.g. `-l 100k`, `-l 1m`.
- `-I` — apply `-l` to inversions (`ori` `><`/`<>`) too; by default they pass
  regardless of distance.
- `-i INT` — maximum `ins_nearby` reads, summed over its sources (default 0,
  i.e. drop any call carrying `ins_nearby`; `-1` = no limit). A call without the
  tag has 0 — `merge` writes it only on calls with `split`/`trim` reads, so calls
  made only of CIGAR `del`/`ins` reads always pass.
- `--src LIST` — comma-separated sources tested by `-c`/`-s` (default: every
  source ending with `.flt`).
- `--rest LIST` — comma-separated sources tested by `-r` (default: every source
  not in `--src`).

Read counts come from the `count=` INFO tag (`count=S.flt:f,r|S.raw:f,r`); each
group's counts are summed over its sources. The typical input merges the raw
breakends (e.g. `S.raw`) with the `fltbnd -s S.flt` survivors, so by default
`-c`/`-s` test the novel `S.flt` reads and `-r 0` drops any call that also has
reads from another source (`S.raw`, or other samples in a multi-sample merge). A line with no source
in either group (e.g. `merge -m` output, whose `count=` carries no source names)
aborts the run, and a warning is printed if no line has a `--src` source at all
(the `*.flt` default doesn't fit — name the sources with `--src`). Survivors are
printed verbatim.

```sh
badclip extract -s S.raw in.bam > S.clip
badclip fltbnd -s S.flt S.clip S.rb3.paf | badclip merge - | badclip fltcnt - > S.sv
```

## Status

All six subcommands are implemented: `extract` (SAM/BAM/CRAM auto-detected, or
PAF via `--paf`) finds breakends; `getfa` turns its output into FASTA; `fltbnd`
filters breakends against a pangenome; `merge` clusters them into consensus SV
calls (with `-m` for a second-level, cross-sample merge); `fltreg` drops
breakends in BED regions; and `fltcnt` keeps calls with enough raw/filtered
support. See `CLAUDE.md` for internals and the interval/offset
convention.
