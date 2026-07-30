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
cargo run -- extract test/join02.paf
```

## Layout

- `src/main.rs`   — clap CLI, subcommand dispatch. New subcommands slot in here.
- `src/io.rs`     — `open_reader`: stdin/file + transparent gzip (magic-byte peek).
- `src/paf.rs`    — `Hit` struct, `Strand`, `parse_paf_line` (keeps only `tp:A:P`).
- `src/extract.rs`— the `extract` command: grouping, sorting, clip/join emission.
- `tests/extract.rs` — end-to-end tests driving the compiled binary.
- `test/`         — `*.paf` inputs, `*.msv` expected outputs, `minisv.js` reference.

## Interval / offset notation

An interval is a pair of **offsets** `|st,en|`, where an offset sits *between*
two adjacent bases, not on a base (0-based, half-open). In `"A|CGT|AGC"`,
`|1,4|` is `"CGT"`. A hit matches a target interval `tc:|ts,te|` to a query
interval `|qs,qe|` on a strand. Emitted positions are raw 0-based offsets — no
±1 adjustment.

## `extract` — behaviour

Input PAF is assumed grouped by read name. Hits are streamed and buffered per
read, keeping only primary/supplementary alignments (`tp:A:P`). Within a read,
hits are sorted by query start `qs`. There is **no** mapq/length filtering
beyond an optional `-q` (default 0, off); the only length knob is `-c` (clip
threshold, default 100).

### Output (9 columns, TAB-delimited)

```
ctg1  pos1  ori  ctg2  pos2  qname  mapq  strand  aln_len=len1,len2
```

`ori` is two characters, each `>`/`<`, or `.` for a missing mate. The final
INFO column carries `aln_len=len1,len2`: the alignment block length (PAF col 11)
of the hit at `ctg1:pos1` and at `ctg2:pos2`, in **output order** (so the values
swap with the endpoints when a join is flipped). For a clip, `len2 = 0`.

- **Left clip** (`first.qs > c`): `first.ctg fts(first) arrow_l.  .  .  qname mapq strand aln_len=alen,0`
- **Right clip** (`last.qlen - last.qe > c`): `last.ctg fte(last) arrow_r.  .  .  qname mapq strand aln_len=alen,0`
- **Join** (each adjacent pair `y0`→`y1`): endpoint `c0 = (y0.ctg, fte(y0), ori(y0), y0.alen)`
  joins `c1 = (y1.ctg, fts(y1), ori(y1), y1.alen)`. Canonicalize so the smaller
  `(ctg, pos)` is first; if swapped, invert both arrows, swap the two `alen`
  values, and set the `strand` column to `-` (else `+`). `mapq` is `min(y0, y1)`.

Helpers: `fts = strand=='+'? ts : te`, `fte = strand=='+'? te : ts`,
`ori = strand=='+'? '>' : '<'`. Clip arrows: `arrow_l = flip(ori)`,
`arrow_r = ori`.

Emission order per read: left clip, joins (read order), right clip.

### Note on `.msv` fixtures

The `test/join*.msv` files carry an INFO column
(`SVTYPE=BND;...;aln_len=...;source=foo`) produced by `minisv.js` whose
`aln_len` differs from badclip's — minisv uses query-span in read order, badclip
uses PAF col-11 block length in output order. So tests compare against each
`.msv` truncated to its first 8 columns with badclip's own `aln_len=...`
appended (see `tests/extract.rs`).

## Not implemented (deliberately, for now)

- BAM input (planned next; will convert BAM records into the same `Hit` and
  reuse `extract.rs`).
- CIGAR-based indel (INS/DEL) extraction.
- Same-contig SV typing (INS/DEL/DUP/INV) and the trailing INFO column — always
  emits `BND`-style orientation regardless of contig.
