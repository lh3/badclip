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
badclip extract [OPTIONS] [PAF]
```

- `PAF` — input alignments in PAF format. Omit it (or pass `-`) to read from
  stdin. Gzip'd input is detected and decompressed automatically, whether from
  a file or stdin.

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `-c`, `--min-clip <INT>` | `100` | Minimum clip length to report a clip breakend. |
| `-q`, `--min-mapq <INT>` | `0` | Drop hits with mapping quality below this value. |

Examples:

```sh
badclip extract aln.paf
badclip extract aln.paf.gz
minimap2 ... | badclip extract          # stdin
gzip -dc aln.paf.gz | badclip extract -  # explicit stdin
```

### Input assumptions

- Alignments are grouped by read name (as produced directly by minimap2).
- Only primary and supplementary alignments (`tp:A:P`) are used; secondary
  alignments are ignored.
- Within a read, hits are sorted by their start position along the read.

## Output

Tab-delimited, 9 columns:

```
ctg1   pos1   ori   ctg2   pos2   qname   mapq   strand   aln_len=len1,len2
```

`ori` is two characters, each `>` or `<` (or `.` for a missing mate).
Positions are raw 0-based offsets that sit *between* bases. The final INFO
column reports `aln_len=len1,len2`, the alignment length (PAF column 11) of the
hit at `ctg1:pos1` and at `ctg2:pos2`; a clip has no second hit, so `len2 = 0`.

- **Clip** — one side is a real breakend, the other is nothing (`.`):

  ```
  chr1   123569841   >.   .   .   read/ccs   1   -   aln_len=14732,0
  chr1   123555134   <.   .   .   read/ccs   1   -   aln_len=14732,0
  ```

- **Join** — two mapped ends meet. The smaller `(contig, position)` is written
  first; the two arrows show how the sides are oriented:

  ```
  chr1    57375269   >>   chr21   32069271   read/ccs   60   -   aln_len=16879,26284
  chr13   51911798   <>   chr2    76026062   read/ccs   60   +   aln_len=29505,5668
  ```

  `chr1:57375269 >> chr21:32069271` means the **right** side of `chr1:57375269`
  joins the **left** side of `chr21:32069271`. `chr13:51911798 <> chr2:76026062`
  means the **left** side of `chr13:51911798` joins the **left** side of
  `chr2:76026062`.

## Status

PAF input is implemented. BAM input is planned and will reuse the same breakend
logic. See `CLAUDE.md` for internals and the interval/offset convention.
