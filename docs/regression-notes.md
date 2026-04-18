# Regression Notes

Detailed internal regression history, debug notes, and performance
snapshots that used to live in `README.md`.

## M8-tier1 smoke-test results — initial (2026-04-17, pre-fix)

First manual pass of the tier-1 matrix on `data/Mp/chr1` x 10 000 PE
reads. Index built by C++ STAR (`--sjdbGTFfile chr1.gtf --sjdbOverhang
99 --genomeSAindexNbases 12`). Matrix: `runThreadN` in `{1,4,8}` x
`outSAMtype` in `{SAM, BAM Unsorted}` x `quantMode` in `{- , GeneCounts}` x
`twopassMode` in `{None, Basic}` — 24 configurations, 12 SAM + 12 BAM.

### What went green

- No crashes or hangs after fixing one `usize::MAX` sentinel bug in
  `crates/star-align/src/windows.rs`.
- Perfect Rust-internal determinism: all three thread counts (1/4/8)
  produce byte-identical `Aligned.out.sam`, `Aligned.out.bam`
  (sorted-body comparison), `SJ.out.tab`, `ReadsPerGene.out.tab`.
- BAM valid: `samtools view` decodes all 12 BAM outputs; 3244
  records each, consistent with the SAM body.
- `tests/e2e.sh` still green (36 cases at that point).
- Faster + smaller footprint than C++ on that fixture
  (`--runThreadN 8 --outSAMtype SAM`):

  | tool | wall | peak RSS |
  | ---- | ---- | -------- |
  | C++ STAR 2.7.11b | 32.5 s | 1.19 GB |
  | star-rs | 12.5 s | 665 MB |

### What diverged from C++ STAR

On that first real-world fixture the Rust port lost coverage:

| metric | star-rs | C++ STAR | delta |
|--------|---------|----------|-------|
| Uniquely mapped reads | 1207 | 1509 | -20 % |
| Multi-mapped reads | 100 | 117 | -14 % |
| `SJ.out.tab` entries | 485 | 740 | -34 % |
| SAM body lines | 3244 | 4034-4036 | -19..20 % |
| SAM body exact matches | 2388 | 2388 | (856 rs-only, 1648 cpp-only) |

This established that the remaining gap was upstream in the seed /
stitch / scoring path, not in sjdb re-insertion.

## M8.0 — representative-read parity (2026-04-17)

Picked one C++-only read —
`A00838:241:H3FVTDSXY:2:1101:10031:7592`, a 150 bp paired-end pair that
C++ maps as a single uniquely-placed `150M`/`150M` alignment at
`chr1:8849598`/`chr1:8849682` with zero mismatches and no splicing.
Rust was reporting the same pair as `unmapped: other`.

### Root cause — `G_OFFSET` / `LOAD_L` mismatch

STAR reserves padding bytes on both sides of the packed genome so that
seed-extension readers can safely spill off either end of a chromosome.
The C++ reference uses:

| path | offset constant | value |
| ---- | --------------- | ----- |
| `Genome::genomeSequenceAllocate` | `Gout = G1out + 100` | 100 |
| `Genome::genomeLoad` | `L = 200` | 200 |

The Rust port had copied both offsets verbatim and then routed all genome
access through `Genome::g_ptr()` using `G_OFFSET`. On a loaded index this
made `g_ptr()` return a pointer 100 bytes into the left pad instead of the
genome start. Downstream, `compare_seq_to_genome` read through the padding
and collapsed good seeds into junk windows.

### Fix

Standardize Rust on `G_OFFSET = 200` for both generation and loading.

### Verification

1. Representative read alone — byte-exact.
2. Full `tests/e2e.sh` regression suite — still green.
3. Full tier-1 10k PE matrix re-run:

   | metric | star-rs (pre-fix) | star-rs (post-fix) | C++ STAR |
   |--------|-------------------|--------------------|----------|
   | Uniquely mapped | 1207 | 1497 | 1509 |
   | Multi-mapped | 100 | 125 | 117 |
   | Unmapped: other | — | 198 | 198 |
   | Unmapped: too short | — | 8180 | 8176 |
   | `SJ.out.tab` entries | 485 | 633 | 740 |
   | SAM body lines | 3244 | 4052 | 4036 |

This closed the catastrophic gap but still left stitch/extension
preference differences.

## M8.0 — follow-up fixes (2026-04-17)

### Fix 2/3: `alignSJstitchMismatchNmax` sentinel

`alignSJstitchMismatchNmax = -1` must behave like "no limit". Rust was
treating it as a literal signed cap instead of C++'s unsigned sentinel.
Fixed in `crates/star-align/src/stitch.rs`.

### Fix 3/3: `alignSJoverhangMin` default

Rust had defaulted `alignSJoverhangMin` to `3`; C++ uses `5`.
Fixed in `crates/star-params/src/parameters.rs`.

### Fix 4/4: `scoreStitchSJshift` + `outSJfilter*`

Two independent issues:

- `scoreStitchSJshift` default was `200`, C++ uses `1`.
- `outSJfilter*` was implemented as a single-pass filter instead of the
  C++ two-stage pipeline.

Fixed in:

- `crates/star-params/src/parameters.rs`
- `crates/star-sjdb/src/output_sj.rs`

## Ath RNA-seq demo regression (2026-04-17)

Real-world sanity run on `RNAseqPluginDemoData/Ath.partial.genome.fna`
(1 MB, 1 chromosome) + its GTF vs. `SRR6281254` chunk from
`RNAseqPluginDemoData/RNAseqData` (251 524 PE x 151 bp, 150 M total bases).
Genome built by C++ STAR with `--sjdbGTFfile ... --sjdbOverhang 100
--genomeSAindexNbases 9`; both aligners run `--runThreadN 8
--outSAMtype SAM --quantMode GeneCounts`.

### Intermediate state after threshold fixes

| metric (Ath 251 k PE, t = 8) | C++ | Rust | delta |
| ---------------------------- | --- | ---- | ----- |
| alignment lines (`Aligned.out.sam` body) | 451 126 | 450 588 | -538 |
| byte-exact SAM-body rate |  |  | 99.885 % |
| Uniquely mapped reads | 221 506 | 221 571 | +65 |
| Splices total | 190 579 | 196 838 | +6 259 |
| Splices non-canonical | 360 | 1 464 | +1 104 |

The remaining divergence was dominated by stitch/extend preference:
Rust sometimes chose a higher-score splice where C++ preferred a
soft-clipped alignment.

### Second sample sanity check — SRR6281256

Same pipeline on `SRR6281256.sra_{1,2}.sample.fastq.gz`:

| metric | C++ | Rust | delta |
| ------ | --- | ---- | ----- |
| alignment lines | 351 750 | 351 732 | -18 |
| byte-exact SAM-body rate |  |  | 99.880 % |
| Uniquely mapped reads | 172 642 | 172 670 | +28 |
| Splices total | 150 559 | 155 605 | +5 046 |
| Splices non-canonical | 314 | 1 243 | +929 |

Same shape as SRR6281254; no new bug class surfaced.

### No-GTF sanity check

To isolate whether sjdb was hiding anything, the Ath genome was
re-indexed without `--sjdbGTFfile` and both samples were rerun:

| sample | byte-exact (no GTF) | byte-exact (with GTF) |
| ------ | ------------------- | --------------------- |
| SRR6281254 | 99.867 % | 99.885 % |
| SRR6281256 | 99.871 % | 99.880 % |

Removing the GTF did not surface a new bug class.

## Ath post-fix — byte-exact at t=1 (2026-04-17)

After the stitch-shift and `outSJfilter*` fixes, single-threaded Ath runs
matched C++ bit-for-bit on every tracked output:

| sample | `Aligned.out.sam` body | `SJ.out.tab` | `ReadsPerGene.out.tab` |
| ------ | ---------------------- | ------------ | ---------------------- |
| SRR6281254 | OK | OK | OK |
| SRR6281256 | OK | OK | OK |

At `--runThreadN 8`, the alignment records remained content-identical
(`sort | md5sum` matches) and `SJ.out.tab` / `ReadsPerGene.out.tab`
remained byte-exact; only record order in `Aligned.out.sam` still differed
because Rust splits reads evenly across threads while C++ chunks by bytes.

## Performance snapshot (2026-04-18)

Measured on the checked-in root `tests/` fixture:

- `genomeGenerate`: `chr1.fa + chr1.gtf`, `--runThreadN 4`,
  `--genomeSAindexNbases 11`, `--sjdbOverhang 149`.
- `alignReads`: `reads_10k_1.fq + reads_10k_2.fq`, `--runThreadN 4`,
  `--quantMode GeneCounts`.
- For `alignReads`, both binaries used the same C++-built index.
- `alignReads` content check on this setup:
  - `Aligned.out.sam` body: byte-exact
  - `SJ.out.tab`: byte-exact
  - `ReadsPerGene.out.tab`: byte-exact
  - `Log.final.out`: stats lines byte-exact after the `mappedReadsU` /
    `transcriptStats` fix

`alignReads` was run 3 times and averaged; peak RSS is the
representative `/usr/bin/time -v` single-run measurement.

| workload | tool | wall | mapping speed | peak RSS |
| -------- | ---- | ---- | ------------- | -------- |
| `alignReads` | C++ STAR 2.7.11b | 30.44 s avg (3 runs) | 1.20 M reads/hour | 745.0 MiB |
| `alignReads` | `star-rs` | 14.98 s avg (3 runs) | 2.40 M reads/hour | 494.7 MiB |
| `genomeGenerate` | C++ STAR 2.7.11b | 13.10 s | — | 920.6 MiB |
| `genomeGenerate` | `star-rs` | 41.90 s | — | 1128.9 MiB |

So on this fixture the current port is:

- `alignReads`: about 2.03x faster than C++ and 33.6 % lower peak RSS.
- `genomeGenerate`: about 3.20x slower than C++ and 22.6 % higher peak RSS.

## Summary of M8.0 fixes

| # | file | what |
| - | ---- | ---- |
| 1 | `crates/star-genome/src/genome.rs` | `G_OFFSET` 100 -> 200 to match `LOAD_L`. |
| 2 | `crates/star-align/src/stitch.rs` | `alignSJstitchMismatchNmax = -1` honored as no-limit. |
| 3 | `crates/star-params/src/parameters.rs` | `alignSJoverhangMin` default 3 -> 5. |
| 4 | `crates/star-align/src/mapped_filter.rs` | Integer truncation for `frac*(Lread-1)` thresholds. |
| 5 | `crates/star-align/src/stitch_window.rs` | Integer truncation for `alignSplicedMateMapLminOverLmate`. |
| 6 | `crates/star-align/src/one_read.rs` | `outFilterMismatchNoverReadLmax` bound fixed. |
| 7 | `crates/star-params/src/parameters.rs` | `scoreStitchSJshift` default 200 -> 1. |
| 8 | `crates/star-sjdb/src/output_sj.rs` | Two-stage `outSJfilter*` pipeline. |
| 9 | `crates/star-sjdb/src/output_sj.rs` | `outSJfilterIntronMaxVsReadN` default fixed. |
