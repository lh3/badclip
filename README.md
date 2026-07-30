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

- `PAF` — input alignments in PAF format. Pass `-` to read from stdin. Gzip'd
  input is detected and decompressed automatically, whether from a file or
  stdin. Running `badclip extract` with no input prints this help instead of
  waiting on stdin.

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `-c`, `--min-clip <INT>` | `100` | Minimum clip length to report a clip breakend. |
| `-q`, `--min-mapq <INT>` | `0` | Drop hits with mapping quality below this value. |
| `-a`, `--min-aln-len <INT>` | `0` | Drop hits whose alignment block length (PAF col 11) is below this value. |

Examples:

```sh
badclip extract aln.paf
badclip extract aln.paf.gz
minimap2 ... | badclip extract -         # stdin
gzip -dc aln.paf.gz | badclip extract -  # stdin (also auto-detects gzip)
```

### Input assumptions

- Alignments are grouped by read name (as produced directly by minimap2).
- Only primary and supplementary alignments (`tp:A:P`) are used; secondary
  alignments are ignored.
- Within a read, hits are sorted by their start position along the read.

## Output

Tab-delimited, 9 columns:

```
ctg1   pos1   ori   ctg2   pos2   qname   mapq   strand   aln_len=..;qlen=..;mapq=..
```

`ori` is two characters, each `>` or `<` (or `.` for a missing mate).
Positions are raw 0-based offsets that sit *between* bases. The final INFO
column carries three tags, all in output order (index 1 ↔ `ctg1:pos1`,
index 2 ↔ `ctg2:pos2`):

- `aln_len=len1,len2` — alignment length (PAF column 11) of each hit; a clip has
  no second hit, so `len2 = 0`.
- `qlen=left,middle,right` — query lengths that sum to the read length. For a
  clip, `right` is the clipped length and `middle` is `0`; for a join, `middle`
  is the (possibly negative) query gap between the two hits.
- `mapq=mapq1,mapq2` — mapq of each hit; for a clip `mapq2 = 0`. (The `mapq`
  *column* above is the smaller of a join's two mapqs; this tag keeps both.)

- **Clip** — one side is a real breakend, the other is nothing (`.`):

  ```
  chr1   123569841   >.   .   .   read/ccs   1   -   aln_len=14732,0;qlen=14925,0,141;mapq=1,0
  chr1   123555134   <.   .   .   read/ccs   1   -   aln_len=14732,0;qlen=14774,0,292;mapq=1,0
  ```

- **Join** — two mapped ends meet. The smaller `(contig, position)` is written
  first; the two arrows show how the sides are oriented:

  ```
  chr1    57375269   >>   chr21   32069271   read/ccs   60   -   aln_len=16879,26284;qlen=16808,1,26163;mapq=60,60
  chr13   51911798   <>   chr2    76026062   read/ccs   60   +   aln_len=29505,5668;qlen=29436,1,5661;mapq=60,60
  ```

  `chr1:57375269 >> chr21:32069271` means the **right** side of `chr1:57375269`
  joins the **left** side of `chr21:32069271`. `chr13:51911798 <> chr2:76026062`
  means the **left** side of `chr13:51911798` joins the **left** side of
  `chr2:76026062`.

## Status

PAF input is implemented. BAM input is planned and will reuse the same breakend
logic. See `CLAUDE.md` for internals and the interval/offset convention.
