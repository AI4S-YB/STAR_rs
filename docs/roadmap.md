# Roadmap

This document tracks the remaining feature work after the current
byte-exact core port and the internal M8 regression plan.

## M7-parked TODO

M8 (large-scale regression) requires feature parity on the BAM/wig/quant
output surface. The following items were deliberately parked at the end
of M7 and must be closed before M8 starts.

### T1. `--outSAMtype BAM Unsorted` byte-exact

- Current: valid BGZF BAM produced via `noodles-bam`; `samtools view` body
  matches C++ after sort, but the BAM file is **not** byte-identical to
  C++ STAR (`alignBAM` uses a hand-rolled BAM record encoder + samtools
  BGZF compressor).
- Work:
  - Port `BAMfunctions.cpp` (record encoding, `packBAM1`, `bam_write1`).
  - Port BGZF block layout used by `BAMoutput` (block size, compression
    level, `bgzf_write` semantics via `libdeflate`).
  - Match record field order (name -> flag -> ref -> pos -> mapq -> cigar ->
    rnext -> pnext -> tlen -> seq -> qual -> tags) and tag ordering
    (`NH HI AS nM MD jM jI ...`).
  - Acceptance: `cmp` on `Aligned.out.bam` for `tiny/` SE + PE + 4-thread.

### T2. `--outSAMtype BAM SortedByCoordinate`

- Current: falls back to Unsorted if requested.
- Strategy: two phases — **algorithm parity first, byte-exact parity
  second**. Phase 2a ships a sort pipeline that clones STAR's algorithm
  1:1 but writes the final BAM through `noodles-bam`; phase 2b swaps
  the writer once T1 lands and thereby earns byte-exact output.

**Phase 2a — algorithm clone + semantic match (this cycle)**

- Port the STAR sort pipeline 1:1 at the algorithm level:
  - `BAMoutput::coordOneAlign` — per-record bin decision, unmapped to
    the last bin, attach `iRead`, maintain `binTotalN` / `binTotalBytes`.
  - `BAMoutput::coordBins` — first-batch boundary estimation
    (`outBAMsortingBinStart`), then redistribute cached records into
    real bins.
  - `BAMbinSortByCoordinate` — sort key
    `(coordinate, iRead, record_offset)`.
  - `BAMbinSortUnmapped` — multi-way merge of per-thread unmapped bin
    files by `iRead`.
  - `bamSortByCoordinate` top-level wiring: bin-size aggregation,
    `--limitBAMsortRAM` check, per-bin parallelism, final bin-order
    concatenation.
- I/O layer: `noodles-bam` writer (not STAR's hand-rolled encoder + BGZF).
  Temp per-bin files use a Rust-native layout
  (`ref_id / pos / i_read / record_bytes`) chosen for debug ergonomics;
  format is fixed for the whole phase, no mid-phase swap.
- `iRead` semantics: stable-within-thread only. Global-ordinal alignment
  with C++ STAR is deferred to phase 2c.
- Phase-2a input source is the existing SAM intermediate, avoiding any
  coupling to T1's native BAM encoder work.
- Honor `--outBAMsortingThreadN`, `--outBAMsortingBinsN`,
  `--limitBAMsortRAM`, plus the `ReadAlign_stitchPieces` /
  `BAMoutput_coordUnmapped` chunk wiring.
- Module layout in `crates/star-bam/`:
  - `bam_output.rs` — `coord_one_align`, `coord_bins`, `coord_flush`.
  - `bam_sort_bin.rs` — `sort_mapped_bin`.
  - `bam_sort_unmapped.rs` — `merge_unmapped_bin`.
  - `bam_sort_coord.rs` — `sort_bam_by_coordinate`,
    `collect_bin_stats`, `check_sort_ram_limit`, `merge_sorted_bins`.
  - `bam_sort_record.rs` — `BinRecordMeta { coord, i_read, file_offset }`.
- Acceptance (phase 2a): compare `samtools view
  Aligned.sortedByCoord.out.bam` against C++ STAR, split by flag 0x4:
  - Mapped (`-F 4`): within each `(ref, pos)` bucket, sort by `qname`
    on both sides and `diff`. Tests the bin/sort algorithm without
    coupling to `iRead`.
  - Unmapped (`-f 4`): `sort | diff` as a multiset. Verifies
    completeness without coupling to `iRead` ordering.
- Unit tests target the algorithm functions directly with hand-built
  records:
  - `coord_bins_matches_cpp_fixture`
  - `mapped_bin_sort_key_matches_cpp`
  - `unmapped_merge_matches_cpp_read_order`

**Phase 2b — byte-exact (gated on T1)**

- Once T1 ships a byte-exact `BAM Unsorted` path (ported `alignBAM` +
  BGZF), swap the phase-2a `noodles-bam` writer for the ported encoder.
  The sort algorithm is unchanged; only the leaf writer moves.
- Acceptance upgrades to `cmp` on `Aligned.sortedByCoord.out.bam`
  (+ `.bai` when `--outBAMindex` is set). This is the acceptance that
  feeds M8 exit criteria #1–#3.

**Phase 2c — optional `iRead` strict alignment**

- If byte-identical unmapped ordering is desired before phase 2b (e.g.
  for tighter semantic diffs during debugging), port STAR's `iRead`
  assignment timing — a global ordinal issued inside `alignBAM` before
  `coordOneAlign` — into the Rust path. Gated on need; not a phase-2a
  blocker.

### T3. `--quantMode TranscriptomeSAM` (+ BAM variant)

- Current: `--quantMode GeneCounts` is byte-exact; `TranscriptomeSAM`
  is not wired.
- Work:
  - Port `Transcriptome::quantAlign` (genome -> transcript coordinate
    transform) and `ReadAlign_outputTranscriptSAM`.
  - Hook into the align pipeline via the existing `quant_hook` closure
    (no new crate boundary needed).
  - Support `--quantTranscriptomeBAMcompression` /
    `--quantTranscriptomeSAMoutput` (already parsed in `star-params`,
    plumbing is stubbed).
  - Acceptance: byte-exact `Aligned.toTranscriptome.out.{sam,bam}` for
    `gctest/`.

### T4. `--outWigType` (wiggle / bedGraph)

- Current: not implemented.
- Work:
  - Port `signalFromBAM.cpp` -> per-strand per-read coverage counters.
  - Honor `--outWigStrand`, `--outWigReferencesPrefix`,
    `--outWigNorm`.
  - Acceptance: byte-exact `Signal.*.bg` / `Signal.*.wig`.

### T5. `--bamRemoveDuplicatesType` (`UniqueIdentical` + `UniqueIdenticalNotMulti`)

- Current: not implemented.
- Work:
  - Port `bamRemoveDuplicates.cpp` (post-sort pass over the BAM chunks).
  - Requires T2 (coord-sorted BAM) as prerequisite.
  - Acceptance: byte-exact `Processed.out.bam`.

### T6. `--quantMode GeneCounts` x strand / spliced-bias counters polish

- Current: 3 counters per gene (unstranded, forward, reverse) match C++.
- Nice-to-have (not a blocker): wire `Log.final.out` lines 30-33
  (`% of reads ...`) derived from Quantifications. Defer unless M8 flags it.

### Exit criteria for M7-parked -> M8

1. T1, T2 (phase 2b), T3, T4, T5 all have at least one byte-exact
   fixture in `tests/e2e.sh`. T2 phase 2a has a separate
   semantic-match fixture (mapped-`diff` + unmapped-multiset) that is
   kept green in parallel until phase 2b supersedes it.
2. Full `tests/e2e.sh` green (current 33 cases + new ones).
3. No use of `std::io::sink()` for mapping output when any BAM/wig flag
   is set (i.e. all output paths are real files, verified by test).

## M8: internal regression on Mp dataset

M8 is an in-house end-to-end regression against a real *Marchantia
polymorpha* (Mp) RNA-seq run hosted on this machine. No ENCODE / public
fixtures are used.

### Dataset

Pinned to `/home/xzg/project/STAR_rs/data/Mp` (not committed to the
repo — too large). Contents:

```text
data/Mp/
├── MpTak_v7.1.fa                # full genome (~240 MB)
├── MpTak_v7.1.gtf               # full annotation (~26 MB)
├── MpTak_v7.1.gff               # GFF3 variant (unused by STAR)
├── MP-S1_1.fq.gz                # PE reads R1 (~2 GB gzipped)
├── MP-S1_2.fq.gz                # PE reads R2 (~2 GB gzipped)
├── chr1/
│   ├── chr1.fa                  # chr1 subset (~31 MB)
│   ├── chr1.gtf                 # chr1 annotation (~5 MB)
│   ├── reads_10k_1.fq           # 10k PE subset, uncompressed
│   └── reads_10k_2.fq
└── chrU/
    ├── chrU.fa                  # chrU subset (~5 MB)
    └── chrU.gtf
```

### Tiers

M8 runs in three tiers. Each tier must be green before the next runs.

### M8-tier1 — chr1 subset x 10k PE reads

- Genome: `data/Mp/chr1/chr1.fa` + `chr1.gtf` (`--sjdbOverhang 99`).
- Reads: `reads_10k_{1,2}.fq` (10k pairs, already plain FASTQ).
- Matrix:
  - `runThreadN` in `{1, 4, 8}`
  - `--outSAMtype` in `{SAM, BAM Unsorted, BAM SortedByCoordinate}`
  - `--quantMode` in `{- , GeneCounts, TranscriptomeSAM, "GeneCounts TranscriptomeSAM"}`
  - `--twopassMode` in `{None, Basic}`
  - `--outWigType` in `{None, bedGraph}`
- Outputs compared: `Aligned.out.sam` body, `Aligned.out.bam`,
  `Aligned.sortedByCoord.out.bam`, `Aligned.toTranscriptome.out.bam`,
  `SJ.out.tab`, `ReadsPerGene.out.tab`, `Signal.*.bg`,
  `Chimeric.out.junction`, `Log.final.out` (stats only, no timing).
- Target runtime: < 60 s per configuration on 8 threads.

### M8-tier2 — full Mp genome x 10k PE reads

- Genome: `data/Mp/MpTak_v7.1.fa` + `MpTak_v7.1.gtf`
  (`--sjdbOverhang 99`).
- Reads: same `reads_10k_{1,2}.fq`.
- Matrix: subset of tier 1 (1-thread and 8-thread, SAM + BAM Unsorted,
  `--quantMode GeneCounts`, `--twopassMode Basic`).
- Purpose: exercise sjdb loading + full genome index with real
  multi-chromosome data.

### M8-tier3 — full Mp genome x full MP-S1 PE library

- Reads: `MP-S1_{1,2}.fq.gz` (tens of millions of pairs). Rust STAR
  can read these directly through extension-based compressed text input
  support; no `zcat` staging is required.
- Matrix: one configuration — `runThreadN=16`,
  `--outSAMtype BAM SortedByCoordinate`,
  `--quantMode "GeneCounts TranscriptomeSAM"`,
  `--twopassMode Basic`, `--outWigType bedGraph`.
- Acceptance:
  - `Aligned.sortedByCoord.out.bam` byte-exact vs C++.
  - `ReadsPerGene.out.tab`, `SJ.out.tab`, `Signal.*.bg` byte-exact.
  - Wall-clock within 1.25x of the C++ reference.
  - Peak RSS within 1.25x of the C++ reference.

### Extra prerequisites beyond M7-parked

**T6. compressed text input / `--readFilesCommand` parity**

- Status: done for native text decompression. STAR-rs now opens plain
  text plus `.gz`, `.xz`, `.bz`, and `.bz2` inputs for
  `--readFilesIn`, `--genomeFastaFiles`, `--sjdbGTFfile`, and
  `--sjdbFileChrStartEnd`.
- Validation: `crates/star-cli/tests/compressed_inputs.rs` compares
  compressed FASTA/GTF/FASTQ inputs against the same plain-text inputs.
- Remaining STAR-compatibility work: port
  `Parameters::readFilesCommandString` -> spawn child process with
  stdout piped into the reader. This is still needed for arbitrary
  commands such as `samtools view`, but no longer blocks M8-tier3
  gzipped FASTQ input.

**T7. `--sjdbGTFfile` in `--runMode genomeGenerate`**

- Status: done. Rust `genomeGenerate` now parses `--sjdbGTFfile`,
  emits `geneInfo.tab` / `transcriptInfo.tab` / `exonInfo.tab` /
  `exonGeTrInfo.tab` / `sjdbList.fromGTF.out.tab`, and writes a
  byte-exact sjdb-enhanced `Genome` / `SA` / `SAindex` on the test
  fixtures used so far.
- Validation:
  - `star-rs/tests/fixtures/gctest` — index files and `GeneCounts`
    outputs are byte-exact vs C++ STAR.
  - root `tests/chr1.fa + chr1.gtf + reads_10k_{1,2}.fq` — index files
    and `ReadsPerGene.out.tab` are byte-exact vs C++ STAR.

### M8 driver layout

```text
tests/
├── e2e.sh
└── m8/
    ├── run.sh
    ├── lib.sh
    ├── tier1_chr1.sh
    ├── tier2_fullmp.sh
    └── tier3_mp_s1.sh
```

Environment knobs:

- `MP_DATA=/home/xzg/project/STAR_rs/data/Mp`
- `STAR_REF=/home/xzg/project/STAR_rs2/STAR/bin/Linux_x86_64/STAR`
- `STAR_RS=target/release/star`
- `MP_THREADS=8`
- `KEEP=1`

### Exit criteria for M8

1. All tier-1 combinations byte-exact vs C++ STAR (SAM/BAM bodies,
   SJ.tab, GeneCounts, TranscriptomeSAM BAM, Chimeric, Signal.bg).
2. Tier-2 matrix byte-exact on the full Mp genome.
3. Tier-3 byte-exact for the full MP-S1 run, with runtime and RSS within
   1.25x of the C++ reference.
4. `tests/m8/run.sh tier1` runs green in under 10 minutes on this host.

### Remaining M8 work

1. M8-tier2 — full Mp genome x 10k PE reads, single SAM config (no
   BAM / GC / 2-pass) as a genome-scale sanity check.
2. M8-tier3 — full Mp genome x full MP-S1 PE library, direct gzip-read
   path. Unblocked by native compressed-input support; pending full run.
