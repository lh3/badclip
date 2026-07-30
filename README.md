# badclip

Extract breakends and structural-variant signals from long-read alignments.

`badclip` reads read-to-reference alignments and reports **breakends**: the
points where a read is soft-clipped away from an alignment (*clips*) and the
junctions between chimeric alignments of the same read (*joins*). Joins reveal
structural rearrangements (translocations, inversions, large indels); clips
flag reads that end abruptly against nothing.

## Install

```sh
cargo build --release
# binary at target/release/badclip
```

## Usage

```sh
badclip extract [OPTIONS] [INPUT]
```

- `INPUT` — alignments in **BAM** by default, or **PAF** with `--paf`. Pass `-`
  to read from stdin. PAF input may be gzip'd (auto-detected). Running `badclip
  extract` with no input prints this help instead of waiting on stdin.

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

## Output

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

Convert `extract` output into a FASTA of the extracted breakend sequences. For
each record that has an `eseq` tag, it writes one FASTA entry named
`readName_idx_leftFlank_rightFlank` (from the `idx` and `elen` tags) with the
`eseq` bases; records without `eseq` are skipped.

```sh
badclip extract aln.bam | badclip geteseq > eseq.fa
badclip geteseq extract.txt            # or a file ("-" / omitted = stdin)
```

Example:

```
>m84039_.../234884533/ccs_0_250_250
ACTTTGGGAGGCCAAGGCAGGCGGATCACCTG...
```

## Status

BAM (default) and PAF (`--paf`) input are implemented and share the same
breakend logic; `geteseq` turns extract output into FASTA. See `CLAUDE.md` for
internals and the interval/offset convention.
