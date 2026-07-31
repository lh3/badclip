> [!Warning]
> This project is vibe coded with Claude Code.

# badclip

Extract breakends and structural-variant signals from long-read alignments.

`badclip` reads read-to-reference alignments and reports **breakends**: the
points where a read is soft-clipped away from an alignment (*clips*) and the
junctions between chimeric alignments of the same read (*joins*). Joins reveal
structural rearrangements (translocations, inversions, large indels); clips
flag reads that end abruptly against nothing.

It has four subcommands, forming a pipeline:

| Subcommand | Purpose |
|------------|---------|
| [`extract`](#extract) | Find breakends in read alignments (the main step). |
| [`geteseq`](#geteseq) | Turn `extract` output into a FASTA of breakend sequences. |
| [`flteseq`](#flteseq) | Drop breakends already explained by a pangenome. |
| [`merge`](#merge) | Cluster per-read breakends into consensus SV calls. |

Every subcommand reads its input from a file or from stdin via `-`, transparently
decompresses gzip'd input, and prints its own `--help` (instead of waiting on
stdin) when run with no input.

## Install

```sh
cargo build --release
# binary at target/release/badclip
```

## `extract`

Find breakends (clips and joins) in read-to-reference alignments.

```sh
badclip extract [OPTIONS] [INPUT]
```

- `INPUT` — alignments in **BAM** by default, or **PAF** with `--paf`. Pass `-`
  to read from stdin. PAF input may be gzip'd (auto-detected).

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `--paf` | off | Read PAF (optionally gzip'd) instead of BAM. |
| `-c`, `--min-clip <INT>` | `50` | Minimum clip length to report a clip breakend. |
| `-q`, `--min-mapq <INT>` | `0` | Drop hits with mapping quality below this value. |
| `-a`, `--min-aln-len <INT>` | `0` | Drop hits whose alignment block length (PAF col 11) is below this value. |
| `-f`, `--flank <INT>` | `250` | Flanking read sequence extracted on each side of a breakend (BAM only). |
| `-e`, `--max-eseq <INT>` | `1000` | Maximum extracted window; longer windows omit `eseq` (BAM only). |

Examples:

```sh
badclip extract aln.bam                   # BAM (default), sorted or unsorted
samtools view -b … | badclip extract -    # BAM from stdin
badclip extract --paf aln.paf.gz          # PAF (gzip auto-detected)
minimap2 … | badclip extract --paf -      # PAF from stdin
```

### Input assumptions

- **BAM**: a read's chimeric hits are read from the primary alignment's `SA:Z:`
  tag; secondary/supplementary/unmapped records are ignored. Works on both
  coordinate-sorted and name-grouped BAM (no grouping required).
- **PAF**: alignments are grouped by read name (as produced directly by
  minimap2); primary/supplementary (`tp:A:P`) and inversion (`tp:A:I`)
  alignments are used, secondary (`tp:A:S`) are ignored.
- Within a read, hits are sorted by their start position along the read.

BAM and PAF produce the same output for the same alignments, with one caveat:
`aln_len` for supplementary hits is a few bp smaller from BAM because minimap2's
`SA` CIGAR is collapsed (the primary hit's `aln_len` matches).

### Output

Tab-delimited, 9 columns:

```
ctg1   pos1   ori   ctg2   pos2   qname   mapq   strand   idx=..;aln_len=..;qlen=..;mapq=..[;elen=..;eseq=..]
```

`ori` is two characters, each `>` or `<` (or `.` for a missing mate).
Positions are raw 0-based offsets that sit *between* bases. The INFO column
begins with `idx`, then three tags in output order (index 1 ↔ `ctg1:pos1`,
index 2 ↔ `ctg2:pos2`):

- `idx=N` — 0-based index of this clip/breakend within the read, in emission
  order (left clip, joins, right clip); resets per read.
- `aln_len=len1,len2` — alignment length (PAF column 11) of each hit; a clip has
  no second hit, so `len2 = 0`.
- `qlen=left,middle,right` — query lengths that sum to the read length. For a
  clip, `right` is the clipped length and `middle` is `0`; for a join, `middle`
  is the (possibly negative) query gap between the two hits.
- `mapq=mapq1,mapq2` — mapq of each hit; for a clip `mapq2 = 0`. (The `mapq`
  *column* above is the smaller of a join's two mapqs; this tag keeps both.)

The last two tags are **BAM only** (PAF has no read sequence), in read-forward
orientation (not flipped with the output):

- `elen=leftFlank,qdist,rightFlank` — sizes of the extracted window: up to `-f`
  bases either side of the junction (`qdist` is the query gap, `0` for a clip),
  so `|eseq| = leftFlank + |qdist| + rightFlank`.
- `eseq=<bases>` — the read sequence over that window, on the **original read
  strand**. Omitted (while `elen` stays) when the window exceeds `-e`.

- **Clip** — one side is a real breakend, the other is nothing (`.`):

  ```
  chr1   123569841   >.   .   .   read/ccs   1   -   idx=0;aln_len=14732,0;qlen=14925,0,141;mapq=1,0
  chr1   123555134   <.   .   .   read/ccs   1   -   idx=1;aln_len=14732,0;qlen=14774,0,292;mapq=1,0
  ```

- **Join** — two mapped ends meet. The smaller `(contig, position)` is written
  first; the two arrows show how the sides are oriented:

  ```
  chr1    57375269   >>   chr21   32069271   read/ccs   60   -   idx=0;aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60
  chr13   51911798   <>   chr2    76026062   read/ccs   60   +   idx=0;aln_len=29505,5668;qlen=29436,1,5661;mapq=60,60
  ```

  `chr1:57375269 >> chr21:32069271` means the **right** side of `chr1:57375269`
  joins the **left** side of `chr21:32069271`. `chr13:51911798 <> chr2:76026062`
  means the **left** side of `chr13:51911798` joins the **left** side of
  `chr2:76026062`.

## `geteseq`

Convert `extract` output into a FASTA of the extracted breakend sequences.

```sh
badclip geteseq [INPUT]
```

- `INPUT` — `extract` output (gzip ok; `-` or omit for stdin).

For each record that has an `eseq` tag, it writes one FASTA entry named
`readName_idx_leftFlank_rightFlank` (from the `idx` and `elen` tags) with the
`eseq` bases; records without `eseq` are skipped.

```sh
badclip extract aln.bam | badclip geteseq - > eseq.fa   # "-" = stdin
badclip geteseq extract.txt                             # or a file
```

Example output:

```
>m84039_.../234884533/ccs_0_250_250
ACTTTGGGAGGCCAAGGCAGGCGGATCACCTG...
```

## `flteseq`

Keep only breakends whose junction sequence is **not** already explained by a
pangenome. Align the `geteseq` FASTA against a pangenome index (e.g.
`ropebwt3 sw`) to get a PAF, then:

```sh
badclip flteseq [OPTIONS] <extractOut> <ropebwt3.paf>
```

- `<extractOut>` — `extract` output (gzip ok; `-` for stdin).
- `<ropebwt3.paf>` — the `ropebwt3 sw` PAF of the `geteseq` FASTA (gzip ok).

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `-l`, `--margin <INT>` | `50` | Margin extending the protected junction interval on each side. |

A record is **dropped** if it has no `eseq` tag, or if some ropebwt3 alignment's
query interval `[qs,qe]` on the eseq contains the junction interval
`[max(0, elen0-l), min(elen0+|elen1|+l, |eseq|)]` (i.e. the alignment spans the
breakend, with `-l` margin on each side — the breakend is "protected"). Surviving
`extract` lines are printed verbatim.

The PAF's query names must be the `geteseq` names, and in the same order as the
`eseq` records (as produced by the pipeline), so the two files are streamed
without loading into memory.

```sh
badclip extract aln.bam > clip.txt
badclip geteseq clip.txt | ropebwt3 sw -d pangenome.fmd - > rb3.paf
badclip flteseq clip.txt rb3.paf > novel.txt
```

## `merge`

Cluster per-read breakends from `extract` into consensus SV calls, each with a
supporting-read count. Unlike the streaming subcommands above, `merge` loads its
input into memory and sorts it by `(ctg1, pos1)`, so the input need not be
pre-sorted.

```sh
badclip merge [OPTIONS] [INPUT]
```

- `INPUT` — `extract` output (gzip ok; `-` or omit for stdin).

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `-c`, `--min-cnt <INT>` | `4` | Minimum supporting-read count to emit a call. |
| `-s`, `--min-cnt-strand <INT>` | `2` | Minimum read count on each strand. |
| `-w`, `--win-size <INT>` | `100` | Clustering window (bp). |
| `-A`, `--max-allele <INT>` | `100` | Cap on simultaneously-open clusters. |
| `-C`, `--max-check <INT>` | `500` | Maximum reads compared per cluster. |

Breakends within `-w` bp on both endpoints (with `><`/`<>` treated as the same
inversion) are grouped; a group is reported only if it has at least `-c` reads
with at least `-s` on each strand. Output keeps the 9-column, TAB-delimited
`extract` shape, with the supporting-read count in the `mapq` slot and a
merge-derived INFO column:

```
ctg1   pos1   ori   ctg2   pos2   .   count   strand   avg_mapq=..;count=..[;count_fr=..;count_rf=..[;foldback]];reads=..
```

- `count=nFwd,nRev` — supporting reads on each strand (the column-7 `count` is
  their sum).
- `count_fr`, `count_rf` — reads oriented `><` and `<>`; present only when the
  cluster has at least one such read. `foldback` is added when only one of the
  two orientations is present (`count_fr` or `count_rf` is `0`) and both
  endpoints are on the same contig.
- `reads=name,…` — the supporting read names.

```sh
badclip extract aln.bam | badclip merge - > sv.txt
```

## Status

All four subcommands are implemented: `extract` (BAM default, or PAF via
`--paf`) finds breakends; `geteseq` turns its output into FASTA; `flteseq`
filters breakends against a pangenome; and `merge` clusters them into consensus
SV calls. See `CLAUDE.md` for internals and the interval/offset convention.
