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
- `src/io.rs`     — `open_raw` (file, or `-` = locked stdin, claimed once per run: a second `-` is an error) under `open_reader` (+ transparent gzip, magic-byte peek; every text input) and `aln::run` (SAM/BAM/CRAM); `line_error`/`display_path` (`path:line: msg` errors, `-` shown as `stdin`).
- `src/record.rs` — shared `extract`-record helpers: `Ctype` (`from_info` requires a known `ctype=`; used by `getfa`/`fltbnd`/`merge`).
- `src/paf.rs`    — `Hit`/`Gap` structs, `Strand`, `parse_paf_line` (keeps `tp:A:P`/`I`; gaps from `cg:Z`), `parse_cigar`, `find_gaps` (long CIGAR indels + insertion normalization).
- `src/aln.rs`    — SAM/BAM/CRAM input via noodles-util (format autodetected): build `Hit`s from a primary record + its `SA` tag; a supplementary record yields only its own gaps. CRAM needs `-r`.
- `src/extract.rs`— dispatch (alignment file vs `--paf`), grouping, sorting, trim/split/del/ins emission; shared `emit_read`/`emit_gaps`/`passes_filter`, `SeqView`.
- `src/getfa.rs`  — `getfa` subcommand (alias `geteseq`): `extract` output → FASTA of eseq windows; `parse_windows` is shared with `fltbnd`.
- `src/fltbnd.rs` — `fltbnd` subcommand (alias `flteseq`): filter breakends by pangenome eseq alignment.
- `src/merge.rs`  — `merge` subcommand: cluster per-read breakends into consensus SV calls.
- `src/fltreg.rs` — `fltreg` subcommand: drop `extract`/`merge` breakends that fall in BED regions.
- `src/fltcnt.rs` — `fltcnt` subcommand: keep `merge` calls by raw/filtered source read counts.
- `tests/extract.rs`, `tests/merge.rs`, `tests/fltreg.rs`, `tests/fltcnt.rs` — end-to-end tests driving the compiled binary.
- `test/`         — `*.paf` inputs, `*.msv` (minisv) references, `bam01*`/`cram01*`/`gap01*`/`flt01*`/`merge01*` fixtures/goldens, `fltreg01.bed`, `minisv.js`.

Cross-subcommand conventions: `main` owns the one `BufWriter` on stdout, passes
it to each subcommand's `run(.., out)` and flushes it once at the end (a final
write error — full disk, closed pipe — is reported, not lost in the writer's
drop); at most one input per run can be `-` (`io::open_raw`); length options
take a `k`/`m`/`g` suffix (`main.rs::parse_num`: `extract -c/-g/-a/-f/-e`,
`fltbnd -l`, `merge -w`, `fltreg -l`, `fltcnt -l`); `getfa`/`fltbnd`/`merge`
require `ctype` (absent — pre-0.2 output — or unknown is a `path:line:` error,
even on a line a filter would drop).

## `getfa`

(Alias `geteseq`, the pre-0.2 name: clap `alias`, hidden from `--help`.)
Filters `extract` output (TAB-delimited; `-` = stdin, gzip auto-detected; no
input → print help) to FASTA. For each eseq window of a record it emits
`>readName_idx_L_R_code` (from the `idx` tag and that window's `elen`
`leftFlank`/`rightFlank`) followed by the bases. `code` is the breakend type:
`T` trim, `S` split, `D` del, and for an insertion's two windows `L` (`eseq`,
its first breakpoint along the read) then `R` (`eseq2`). A record without
`ctype` (pre-0.2 output) or with an unknown one is an error (`Ctype::from_info`);
a line with < 9 columns is skipped. Records lacking an eseq are skipped. `-Q` (default 20)
additionally drops records whose `equal` quality — the lower of
`equal`/`equal2`, so an insertion emits both windows or neither — is below the
threshold; records without an `equal` tag are kept. The parser
(`getfa.rs::parse_windows`) is shared with `fltbnd` so names always agree.
Golden: `test/bam01.fa`.

## `fltbnd`

(Alias `flteseq`, hidden like `geteseq`.) Filters `extract` output against a
ropebwt3 `sw` PAF (both gzip-ok) of the `getfa` FASTA aligned to a pangenome.
Drops a line if it has no eseq, or if some PAF alignment's query interval
`[qs,qe]` (cols 3,4) contains the junction interval
`[max(0,e0-l), min(e0+|e1|+l, e0+|e1|+e2)]` on one of its eseq windows
(`elen=e0,e1,e2`, `-l` margin, default 50) — the breakend is "protected"/not
novel. An insertion has two windows (`_L`, `_R`) and is dropped if **either**
junction is protected. Survivors are printed verbatim. The PAF qname is
`getfa`'s name and its lines are in the same order as the FASTA records (a
subset), so both files stream with a one-line PAF lookahead — no in-memory load;
every window's PAF lines are consumed even once one is protected, to stay in
sync. At the end, a PAF line left unconsumed (blank lines skipped,
`fltbnd.rs::next_paf`) is an **error**: it matched no window in order (input
filtered after `getfa`, PAF reordered), so every breakend after the stall was
judged novel. A window `getfa -Q` dropped has no PAF lines and is novel to a
lower `fltbnd -Q` — undetectable, so the docs say to use the same `-Q`. Either
input may be `-` (stdin), not both (`io::open_raw`); missing input → help (via
`main.rs::print_subcommand_help`, shared with `extract`/`getfa`). Fixtures:
`test/flt01.clip`, `test/flt01.rb3.paf`.

`-Q` (default 20) additionally drops any line whose `equal` quality (lower of
`equal`/`equal2`) is below the threshold (in both modes; lines without an
`equal` tag are kept). The check runs **after** the PAF lookahead consumes that
line's entries, so the stream stays in sync regardless of how the PAF was built.

`-s STR` flips the mode: print all input lines (except those `-Q` drops) and
rewrite the `source=` tag of the INFO column (col 9, not the last column) to
`STR` on the survivors (kept/novel), leaving dropped lines (no `eseq`, or
pangenome-explained) verbatim — tagging novel calls distinctly for a downstream
`merge` (`fltbnd.rs::relabel_source`). Without `-s`, only survivors are printed
verbatim.

## `merge`

Collapses per-read `extract` breakends into consensus SV calls — a simplified
`gc_cmd_merge` from `minisv.js`. badclip has only the breakend type (no
`SVTYPE`/`SVLEN`) and no centromere/RT annotation, so
minisv's SVTYPE/SVLEN checks, the `-d` filter, the centromere filter (`-e`), and
the RT branch (`-r`/`-R`) are all dropped. Unlike minisv (which needs an upstream
`sort -k1,1 -k2,2n`), `merge` **loads all records into memory and sorts them
itself** by `(ctg, pos)` — Rust `str` order == `LC_ALL=C sort`; the sort is
stable so the representative pick is deterministic. Takes one or more inputs,
read in turn as if concatenated (each `-`=stdin, gzip auto-detected per file);
no input → help.

Algorithm (`src/merge.rs`): sweep a window of active clusters. Each record joins
the active cluster whose members it most often matches (`same_sv`), else starts a
new one; a cluster is flushed (and emitted if it passes the filters) once the
sweep moves past its `pos_max + -w`. `same_sv` requires: same `ctg`; compatible
`ori` (equal, or the inversion pair `><`↔`<>`); `|pos−pos| ≤ -w`; same `ctg2`;
and `|pos2−pos2| ≤ -w` **skipped when either `pos2` is absent** (a clip's `.`
mate — reproduces minisv's `NaN` comparison, so clips cluster on `ctg2`/`pos`
alone). Per-cluster comparisons are capped at `-M` members (a deterministic cap
replacing minisv's reservoir sampling — deep-cluster counts are capped, not
scaled; minisv's `-C`, which badclip uses for the `-m` per-line filter), and the
active-cluster list is bounded by `-A`. Clustering ignores `ctype` (`Rec::ctype`,
`record::Ctype`; required — `parse_rec` errors with `path:line:` on a missing or
unknown one before any filter): a CIGAR `del` clusters with split-read deletions
at the same coordinates, and the call reports the mix (`ctype_count`).

Flags mirror minisv's kept subset: `-c` min read count (2), `-s` min count on
each strand (0), `-w` window bp (100), `-A` max active clusters (100), `-M` max
reads compared per cluster (500). `-Q` (default 20) drops input breakends whose
`equal` quality (the lower of `equal`/`equal2`) is below the threshold before
clustering (`parse_rec`); breakends without an `equal` tag are kept. `-q` (default 0) drops input breakends whose
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
inversion `><`/`<>` members average correctly; `q2=0` for a clip cluster), then
`;ctype_count=<type:n|...>` (members per ctype, alphabetical, present types
only; sums to col 7), then — when the cluster has `ins` members —
`;avg_inslen=N;diff_inslen=M` (rounded mean and `max-min` of their `inslen`),
then, when any member is `><`/`<>`, `;count_fr=A;count_rf=B` and `;foldback` if
`A*B==0 && ctg1==ctg2`, then `;ins_nearby=<src:n|...>`, and finally
`;reads=<src:name,...|...>`. `count=` lists each
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
`test/merge01.clip` (unsorted input), `test/merge01.expected` (`-c 3 -s 1`
golden).

`ins_nearby` (only for a cluster with a `split`/`trim` member; omitted when 0)
flags spurious translocations from long germline insertions whose sequence maps
elsewhere: the distinct `ins` reads — deduplicated by `(source, qname)` across
both ends, counted per source, alphabetical — of the **whole input** (any
source, after the per-record `-Q`/`-q`/`-p` filters, incl. those whose own
cluster fails `-c`) near either printed end `(ctg1,pos1)`/`(ctg2,pos2)`, **minus
the call's own reads** (any member's `(source, qname)`, whatever its ctype — so a
split+`ins` insertion call doesn't flag itself, and a split read that also
carries a CIGAR `ins` nearby isn't counted). An `ins`
record `ctg L >> ctg R` is near `(ctg,p)` iff `L-w <= p <= R+w` (`-w`).
`InsIndex` (built in `run` after parsing, before the sweep moves records into
clusters) holds a per-contig `IITree` of `[L-w, R+w+1)` per `ins` record; an end
is a point query `overlap(p, p+1)`.

### `-m` sample-merge

`-m` re-runs the **same clustering algorithm** on the (filtered) combination of
per-sample `merge` outputs to produce population-level, sample-merged calls.
`merge` and `extract` output share the 9 columns but differ in INFO, so `-m`
swaps only the input parser (`parse_rec_merge`) and output emitter
(`write_sv_merge`); `same_sv`/sweep are reused verbatim.

- **Input parsing** (`-m`): mapq from `avg_mapq=q1,q2` (col 7 is a count here,
  not mapq); per-strand read counts from `count=` summed over its
  `|`-separated `[src:]f,r` entries (source-prefix optional via `parse_count`, so
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
  **No `reads=` tag** (too many at 3202-sample scale), and no
  `ctype_count`/`avg_inslen`/`diff_inslen`/`ins_nearby` (ignored on input; no
  `InsIndex` is built). The output is itself
  `-m`-parseable, so merges chain (col-7 then falls back to member count, since
  `count=F,R` drops the source labels).

The `-C`/`-S` per-line filters are no-ops in the default (non-`-m`) path; there
`-Q`/`-q`/`-p` test the `equal`/col-7/`mapq=` values described above.

## `fltreg`

Drops `extract`/`merge` output lines whose breakends land in BED regions.
`fltreg [-l NUM] <in> <bed>` (both positional; either may be `-` = stdin, not
both; gzip auto on both; either missing → help). Since both output formats share the 9 columns, one filter
serves both: a line is dropped if **either** endpoint — `(ctg1,pos1)` (cols 0,1)
or, when present, `(ctg2,pos2)` (cols 3,4) — falls in a region; a clip
(`ctg2="."`) tests only endpoint 1. Survivors are printed verbatim; lines whose
endpoints can't be parsed are kept.

Positions are raw 0-based offsets and BED is half-open, so offset `p` is inside
`[start,end)` iff `start <= p < end`. `-l` (default 0) pads every region to
`[start-l, end+l)`; the padding is applied in `Regions::load` **before** the
merge, so `contains` is unchanged. BED loading (`Regions::load`) skips blank /
`#` / `track` / `browser` lines, takes cols 0/1/2 as `chrom start end`, and per
contig **sorts by start and merges overlapping/adjacent intervals** into a
disjoint list; `Regions::contains` then does one `partition_point` binary search
(the rightmost interval with `start <= p` is the only one that can contain `p`).
No new dependency. Fixture: `test/fltreg01.bed`.

## `fltcnt`

`fltcnt [opts] <in>` keeps `merge` (non-`-m`) lines by read support, printing
survivors verbatim. Counts come from the `count=` tag's `src:f,r` entries, split
into two groups: `--src` (comma list; default every source ending with `.flt`)
and `--rest` (comma list; default every source **not** in `--src`); each group's
`f,r` is summed over its sources (a source in both lists counts in both). A line
with no source in either group (incl. source-less `-m` output) **aborts**; if no
line has any `--src` source, a stderr warning is printed at the end. Keep iff
`--src` total `>= -c` (3), `--src` `min(f,r) >= -s` (0), and `--rest` total
`<= -r` (0; negative = no cap) — `-r` is a **maximum**, so by default a call
survives only if all its reads are from `--src` sources. `-l`
(default 0; `k`/`m`/`g` suffix ok via `main.rs::parse_num`): when `>= 0`, drop
clips (`ctg2 = "."`), and when positive also same-contig lines with
`pos2 - pos1 < -l`; inter-contig joins pass. (An insertion's `pos2 - pos1` is
only its normalization range `R - L`, so any positive `-l` drops insertions.) Inversions (`ori` `><`/`<>`) are
exempt from the distance test (always kept, regardless of `pos2 - pos1`) unless
`-I` is given, which applies `-l` to them too; clips are still dropped either
way. `-l -1` (any negative) disables it. `-i` (default 0) also requires the
`ins_nearby=src:n|...` counts, summed over all sources (absent tag = 0), to be
`<= -i` — a **maximum** like `-r` (negative = no cap), so by default any call
carrying `ins_nearby` is dropped (`fltcnt.rs::ins_nearby`). `merge` writes the
tag only on calls with a `split`/`trim` member (and not counting the call's own
reads), so pure `del`/`ins` calls always pass `-i`. Tests: `tests/fltcnt.rs`.

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
  (inversion), dropping `tp:A:S`. One hit per line; its gaps come from `cg:Z`.
- **Alignment file** (`src/aln.rs::run`): the container (SAM/BAM/CRAM) is
  autodetected by `noodles_util::alignment::io::reader` (magic-byte sniff, works
  on stdin too); records come as `Box<dyn sam::alignment::Record>` and are
  processed by the trait (one path for all three). Trims/splits come from
  primary records only (secondary records are skipped, unmapped ones only feed
  `Stats`); supplementary records yield just their own del/ins (see below).
  Each read's hits = the primary (from its own CIGAR) + one per `SA:Z:` entry.
  For paired-end reads the qname gets a
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
    the `SA` tag are dropped before breakends form (a supplementary record on an
    ALT contig, or whose primary — first `SA` entry — is, is skipped). No-op for
    `--paf`. `--alt -` reads stdin, refused when the input is `-` too (both
    would lock stdin and hang).
  - noodles crates: util 0.82 (feature `alignment`) / sam 0.87 / core 0.20 /
    cram 0.96 / fasta 0.64.

Within a read, hits are sorted by query start `qs`. Filtering is off by default.
`-a` drops hits whose alignment block length (PAF col 11, the `alen` field) is
below a threshold (default 0) — the one **upfront** hit filter (`passes_filter`),
applied before breakends are formed. Trims (`ctype=trim`, read ends clipped away
from every alignment) are emitted only with `-T`; their threshold is `-c`
(default 50). `-s` sets the `source=` dataset label (default `foo`), stamped on
every record (alignment file and PAF).

**CIGAR indels (`ctype=del`/`ins`).** Every `I`/`D` op of at least `-g` bases
(default 1000; `k` suffix ok; `<= 0` disables) becomes a breakend
(`paf.rs::find_gaps`, stored as `Hit::gaps`). The `SA` CIGAR is collapsed, so
SA-derived hits carry no gaps; instead each **supplementary record** is read on
its own (secondary records are skipped) and emits only its own gaps
(`aln.rs::emit_record`): it is skipped if its contig is ALT, if the SA's first
entry (its primary, by SAM convention) is ALT, or if its own hit fails `-a`;
its `n_aln` is self + non-ALT SA hits passing `-a`; it never touches `Stats`;
gaps are found before its SEQ is decoded, so gap-less records stay cheap. PAF
gaps come from `cg:Z` (absent → none). With `q` the event's query offset in
alignment orientation (counting the leading clip) and `t` its reference offset,
a gap is written like the split an aligner would report for it — `ori` `>>`,
col 8 = the alignment strand, the one hit standing in for both sides
(`aln_len=A,A`, `mapq=M,M`, col-7 `M`, `-q` against `M`):

- del: `ctg t >> ctg t+len`, `qlen=q,0,qlen-q`; one window at the read-forward
  junction `qf = +? q : qlen-q`.
- ins: `ctg L >> ctg R`, `qlen=q,len,qlen-q-len`, `inslen=len`; windows at the
  read-forward breakpoints `qa = +? q : qlen-q-len` (`elen/equal/eseq`) and
  `qa+len` (`elen2/equal2/eseq2`), both on the CIGAR placement.

`idx` for a gap is its read-forward offset (`qf`/`qa`): a read's gaps come from
several records in any file order, so they can't share the trim/split counter,
but the offset is unique per read and ctype (so `getfa` names are unique).
Emission order: a primary record emits trims/splits, then each hit's gaps in
read-forward order (only its own hit has gaps in BAM; every hit in PAF).

**Insertion normalization** (alignment input only; `find_gaps` with the stored
SEQ in reference-forward orientation and the leading hard-clip offset). With the
inserted bases at SEQ `[a,b)`: shift right `k_r` while `seq[a+i]==seq[b+i]`,
left `k_l` while `seq[a-1-i]==seq[b-1-i]`, each bounded by the run of
consecutive `M/=/X` ops adjacent to the `I` on that side (so each shifted base
is one reference base) and stopping at an `N` (case-insensitive compare) or the
SEQ end. `L = t-k_l`, `R = t+k_r`; a shift may exceed `len` (tandem repeats).
PAF has no sequence, so its insertions stay at `t >> t`. Deletions are never
shifted (that needs the deleted reference bases). Unit tests:
`paf::tests::insertion_*`.

Two **post-filters** drop emitted lines (never alignments) — like their `merge`
namesakes — so the per-read `idx`/`n_aln` counters are unchanged (dropped
breakends are simply absent, so `idx` stays a stable per-breakend identifier;
likewise a trim hidden by the lack of `-T` still takes its `idx`):
`-Q` (default 0 = keep all) drops breakends whose computed `equal` quality is
below the threshold — for an insertion, if either window's is (breakends
without an `equal` value — no `eseq`, no `QUAL`, or PAF — are always kept);
`-q` (default 0) drops lines whose col-7 mapq (the split's `max`, or the one
hit's mapq for trim/del/ins) is below the threshold. Both shrink the output
feeding downstream tools. In the `--help` listing `-q` is placed after `-Q` to
reflect that it is a post-filter.

**End-of-run stats** (`extract.rs::Stats`, printed to stderr after a successful
run): number of reads, number of unmapped reads, total bases in primary
alignments, and read N50. A read is one primary record (alignment input;
secondary/supplementary records are ignored; unmapped records are included —
their `SEQ` length counts toward N50, with no primary bases; mapped reads are
counted after the ALT skip) or one PAF read group (an unmapped PAF read is a line
whose target name, col 6, is `*` — minimap2 `--paf-no-hit`;
`paf.rs::parse_paf_unmapped`); counted **before** the `-a` filter. Primary bases sum the primary
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
`[max(lo-f,off), min(hi+f,off+|SEQ|)]` where `lo=min(y0.qe,y1.qs)`, `hi=max(..)`
for a split (`lo=hi=p` otherwise: the clip point `first.qs`/`last.qe` for a trim,
the junction for a del, each breakpoint for an ins). `-f` (default 250) is the
flank, `-e` (default 1000) the max window (per window).
`elen=leftFlank,qdist,rightFlank` is always emitted; `eseq=<bases>` is appended
only when the window `<= -e`. Both are read-forward (NOT flipped with the output,
unlike `qlen`). The record's `SEQ`/`QUAL` become a `SeqView` (`aln.rs::seq_view`):
reverse-complemented (qualities reversed) iff the record is reverse-mapped
(`src/aln.rs::revcomp`), covering `|off, off+|SEQ||` — `off=0` when
`SEQ.len()==qlen` (soft-clipped), `off=qs` when `SEQ.len()==qe-qs` (hard-clipped,
e.g. a supplementary record), else no view. Trim/split windows need the whole
read (`off==0`), exactly as before; del/ins windows are clamped to the view.
For CRAM the `SEQ` is reconstructed against the `-r` reference. See
`src/extract.rs::eseq_info`.

`equal=<phred>` is the quality of the eseq window, emitted **between `elen` and
`eseq`** — only when `eseq` is emitted and per-base qualities cover the window
(omitted when `QUAL` is `*`). It turns each base quality to an error rate, means
the expected errors over the window, and scales back to phred:
`round(-10*log10((Σ 10^(-Q_i/10))/L))`, `L=|eseq|` (`eseq_qual`). Qualities are
read-forward too: reversed (not complemented) on a reverse-mapped primary. For
constant quality `Q` this is `Q`.

### Alignment file vs PAF parity

For the same alignments the output is identical except (1) `aln_len` of
supplementary (SA-derived) hits, which is a few bp smaller from a BAM/CRAM
because minimap2's `SA` CIGAR merges small indels, and of hits with long indels,
where minimap2's PAF col 11 differs from the CIGAR's M+I+D; and (2) insertion
positions, normalized only for alignment input (`L <= t <= R`, PAF `t >> t`).
Everything else — line count, offsets, `ori`, `strand`, `qlen`, `mapq`, `idx`,
trim/split/del/ins logic, and inversion (`tp:A:I`) segments — matches exactly.
(Verified on HG002.HiFi-10k.mm2: 1365 lines with `-T` from both, incl. 230 del /
73 ins; diffs only as above. `tests/extract.rs::gap01_bam_paf_parity`.) BAM and
CRAM outputs are byte-identical for the same reads (`tests/extract.rs::cram_matches_bam`).

### Output (9 columns, TAB-delimited)

```
ctg1  pos1  ori  ctg2  pos2  qname  mapq  strand  source=..;ctype=..;idx=..;n_aln=..;aln_len=..;qlen=..;mapq=..[;inslen=..][;elen=..[;equal=..;eseq=..]][;elen2=..[;equal2=..;eseq2=..]]
```

`ori` is two characters, each `>`/`<`, or `.` for a missing mate. Tag **order in
the INFO column is not significant** — every consumer parses by key, not
position — but the emission order is `source`, `ctype`, `idx`, `n_aln`, then
three tags keyed to **output order** — index 1 tracks `ctg1:pos1`, index 2
tracks `ctg2:pos2` — so their values swap with the endpoints when a split is
flipped (`source` is per-read and `elen`/`eseq`, appended for BAM, do not swap):

- `source=NAME` — the `-s` dataset label (default `foo`), for telling reads of
  different datasets apart (tumor/normal, trio) downstream; written first purely
  for readability. Same on every record of a run.
- `ctype=TYPE` — `trim` (`-T` only), `split`, `del` or `ins` (see above).
- `idx=N` — for trim/split, 0-based index within the read in emission order
  (left trim, splits, right trim), resetting per read (`emit_read`'s local
  counter; hidden trims still count); for del/ins, the read-forward offset.
- `n_aln=N` — number of alignments (primary + supplementary) for the read that
  passed the `-a` filter, i.e. `hits.len()` at `emit_read`. A read-level
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

- **Left trim** (`-T`, `first.qs > c`): `first.ctg fts(first) arrow_l.  .  .  qname mapq strand INFO`
- **Right trim** (`-T`, `last.qlen - last.qe > c`): `last.ctg fte(last) arrow_r.  .  .  qname mapq strand INFO`
- **Split** (each adjacent pair `y0`→`y1`): endpoint `c0 = (y0.ctg, fte(y0), ori(y0), y0.alen)`
  joins `c1 = (y1.ctg, fts(y1), ori(y1), y1.alen)`. Canonicalize so the smaller
  `(ctg, pos)` is first — on a tie (same ctg and pos) swap iff `ori(y0)` is `<`,
  so a touching junction is `>>`/`+` from a `+` read and `>>`/`-` from a `-`
  read (minisv swaps on every tie, giving `<<` vs `>>`, which never cluster);
  if swapped, invert both arrows, swap the two `alen`
  values, the two `mapq` values, and the `qlen` `left`/`right`, and set the
  `strand` column to `-` (else `+`). The field-7 `mapq` column is `max(y0, y1)`.

Helpers: `fts = strand=='+'? ts : te`, `fte = strand=='+'? te : ts`,
`ori = strand=='+'? '>' : '<'`. Clip arrows: `arrow_l = flip(ori)`,
`arrow_r = ori`.

Emission order per primary record: left trim, splits (read order), right trim,
then del/ins (read order). A supplementary record emits only its del/ins.

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
`test/inv01.paf` is a `tp:A:I` inversion read (two joins). `bam01.expected` is
made with `-T`; `tests/extract.rs::trims_off_by_default` checks the default
output is it minus the trim lines. `test/gap01.bam` (+ `gap01.paf`, golden
`gap01.expected`) holds every record of four HG002 reads with ≥1 kb CIGAR
indels, incl. hard-clipped reverse supplementary records; its positions, shifts,
`qlen` and eseq/equal were checked independently against `samtools view` and the
FASTQ. `tests/extract.rs::gap_sam` builds synthetic del/ins reads with exact
expected lines.

## Not implemented (deliberately, for now)

- Same-contig SV typing (INS/DEL/DUP/INV) and the trailing INFO column — always
  emits `BND`-style orientation regardless of contig.
