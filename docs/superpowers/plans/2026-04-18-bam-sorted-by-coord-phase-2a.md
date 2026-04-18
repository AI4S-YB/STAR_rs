# BAM SortedByCoordinate (Phase 2a) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship `--outSAMtype BAM SortedByCoordinate` end-to-end, cloning STAR's bin/sort/merge algorithm 1:1 at the algorithm layer, while writing the final BAM through `noodles-bam` (not STAR's hand-rolled encoder). Acceptance is semantic equivalence via `samtools view`, split by `flag 0x4` (mapped: sort-within-coord-bucket diff; unmapped: multiset diff).

**Architecture:** Reuse the existing `Aligned.out.sam.tmp` SAM intermediate as the phase-2a record source — no coupling to T1's native BAM encoder. A single-threaded feeder parses each SAM record, converts to `noodles_bam::Record`, assigns a monotonic `iRead`, and calls `BamOutput::coord_one_align` on a single per-pass `BamOutput` (iChunk=0). Per-bin sorting and unmapped multi-way merge run in parallel across bins via `rayon`. Final assembly writes the output BAM via `noodles_bam::io::Writer`, concatenating sorted per-bin temp records in ascending bin order.

**Tech Stack:** Rust 2021, `noodles-bam 0.88`, `noodles-sam 0.84`, `noodles-bgzf 0.46`, `rayon 1.10`, `anyhow`, `thiserror`. Temp-file layout is a Rust-native binary format (see Task 1), not STAR's `bam_bytes + iRead` layout.

**Non-goals (deferred to phase 2b/2c):** byte-exact `Aligned.sortedByCoord.out.bam`, `iRead` global-ordinal alignment with C++ STAR, `--outBAMindex` (`.bai`), per-thread BamOutput feeding (phase 2a uses a single feeder thread).

---

## File Structure

All changes scoped to two crates: `star-bam` (new modules + struct types), `star-cli` (wiring). `star-params` gains three new fields plus three CLI flags.

**Created files:**
- `crates/star-bam/src/bam_sort_record.rs` — `BinRecordMeta` sort-key struct + temp-file record codec (`TMP_RECORD_MAGIC`, `read_tmp_record`, `write_tmp_record`).
- `crates/star-bam/src/bam_output.rs` — `BamOutput` struct, `coord_one_align`, `coord_bins`, `coord_flush` (replacing current 5-line stub).
- `crates/star-bam/src/bam_sort_bin.rs` — `sort_mapped_bin(ibin, bin_n, bin_s, n_threads, tmp_dir) -> Result<()>` (new module).
- `crates/star-bam/src/bam_sort_unmapped.rs` — `merge_unmapped_bin(ibin, n_threads, tmp_dir) -> Result<()>` (new module).
- `crates/star-bam/src/bam_sort_coord.rs` — `sort_bam_by_coordinate`, `collect_bin_stats`, `check_sort_ram_limit`, `assemble_final_bam` (replacing current 5-line stub).
- `crates/star-bam/src/bam_sort_feed.rs` — `feed_sam_into_bam_output(sam_path, bam_output)` converts SAM intermediate records into `BamOutput::coord_one_align` calls (new module).
- `crates/star-bam/tests/bam_output_algorithm.rs` — unit tests against hand-built `noodles_bam::Record`s.

**Modified files:**
- `crates/star-bam/src/lib.rs` — register new modules.
- `crates/star-bam/Cargo.toml` — add `rayon`, `noodles-core` (already transitively pulled, make explicit).
- `crates/star-params/src/parameters.rs` — add `out_bam_coord_nbins: u32` (default 50), `out_bam_sorting_thread_n: u32` (default 0 → auto), `limit_bam_sort_ram: u64` (default 0 → "use all available"), `out_bam_sorting_bin_start: Vec<u64>` (runtime-populated), `out_bam_sort_tmp_dir: String` (derived from `--outTmpDir`), `chunk_out_bam_size_bytes: u64` (default 1 << 24 = 16 MiB). Parse the three CLI flags.
- `crates/star-cli/src/main.rs:259-263, 394-405` — split the `want_bam` branch into Unsorted vs SortedByCoordinate; the latter calls `bam_sort_coord::sort_bam_by_coordinate`.
- `tests/e2e.sh` — add test case "[9] alignReads --outSAMtype BAM SortedByCoordinate" (single-thread + 4-thread), using mapped `-F 4` sort-within-bucket diff + unmapped `-f 4` multiset diff.

---

## Task 1: `BinRecordMeta` + temp-file record codec

**Files:**
- Create: `crates/star-bam/src/bam_sort_record.rs`
- Create: `crates/star-bam/tests/bam_output_algorithm.rs`
- Modify: `crates/star-bam/src/lib.rs:24` (register module)

**Temp-file record layout** (phase 2a — Rust-native, not STAR byte-for-byte):

```
offset  size  field
0       4     magic = 0x52_42_52_31 ("RBR1")
4       4    record_len (u32, little-endian) — length of bam_bytes
8       4    ref_id (u32) — bam::Record::reference_sequence_id, or 0xFFFF_FFFF for unmapped
12      4    pos (u32) — 0-based alignment start, or 0xFFFF_FFFF for unmapped
16      8     i_read (u64, little-endian)
24      record_len   bam_bytes — raw BAM record bytes (including the 4-byte block_size prefix)
```

- [ ] **Step 1: Write failing tests**

Create `crates/star-bam/tests/bam_output_algorithm.rs`:

```rust
//! Unit tests for the Rust STAR BAM sort algorithm (phase 2a).
//! These target algorithm-level parity with STAR's C++ (BAMoutput.cpp +
//! BAMbinSortByCoordinate.cpp + BAMbinSortUnmapped.cpp) using hand-built
//! `noodles_bam::Record`s; they do not depend on a real STAR run.

use star_bam::bam_sort_record::{read_tmp_record, write_tmp_record, BinRecordMeta, TmpRecord, UNMAPPED_SENTINEL};

#[test]
fn tmp_record_roundtrip_mapped() {
    let mut buf: Vec<u8> = Vec::new();
    let bam_bytes = vec![0xAA, 0xBB, 0xCC, 0xDD];
    let tmp = TmpRecord { ref_id: 2, pos: 1000, i_read: 7, bam_bytes: bam_bytes.clone() };
    write_tmp_record(&mut buf, &tmp).unwrap();
    assert_eq!(buf.len(), 24 + bam_bytes.len());
    let mut cursor = std::io::Cursor::new(&buf);
    let got = read_tmp_record(&mut cursor).unwrap().expect("record");
    assert_eq!(got.ref_id, 2);
    assert_eq!(got.pos, 1000);
    assert_eq!(got.i_read, 7);
    assert_eq!(got.bam_bytes, bam_bytes);
}

#[test]
fn tmp_record_roundtrip_unmapped() {
    let mut buf: Vec<u8> = Vec::new();
    let tmp = TmpRecord { ref_id: UNMAPPED_SENTINEL, pos: UNMAPPED_SENTINEL, i_read: 42, bam_bytes: vec![1, 2, 3] };
    write_tmp_record(&mut buf, &tmp).unwrap();
    let mut cursor = std::io::Cursor::new(&buf);
    let got = read_tmp_record(&mut cursor).unwrap().expect("record");
    assert_eq!(got.ref_id, UNMAPPED_SENTINEL);
    assert_eq!(got.pos, UNMAPPED_SENTINEL);
    assert_eq!(got.i_read, 42);
}

#[test]
fn tmp_record_eof_returns_none() {
    let buf: Vec<u8> = Vec::new();
    let mut cursor = std::io::Cursor::new(&buf);
    let got = read_tmp_record(&mut cursor).unwrap();
    assert!(got.is_none());
}

#[test]
fn bin_record_meta_sort_order() {
    // (coord, i_read, offset) sort: primary coord ascending, then iRead ascending.
    let mut metas = vec![
        BinRecordMeta { coord: 100, i_read: 5, file_offset: 1000 },
        BinRecordMeta { coord: 50,  i_read: 9, file_offset: 2000 },
        BinRecordMeta { coord: 100, i_read: 3, file_offset: 3000 },
        BinRecordMeta { coord: 50,  i_read: 1, file_offset: 4000 },
    ];
    metas.sort_by(|a, b| a.sort_key().cmp(&b.sort_key()));
    let coords: Vec<_> = metas.iter().map(|m| (m.coord, m.i_read)).collect();
    assert_eq!(coords, vec![(50, 1), (50, 9), (100, 3), (100, 5)]);
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p star-bam --test bam_output_algorithm 2>&1 | head -40`
Expected: compilation error — module `bam_sort_record` not found.

- [ ] **Step 3: Create `bam_sort_record.rs`**

Create `crates/star-bam/src/bam_sort_record.rs`:

```rust
//! Temp-file record codec for the phase-2a BAM sort pipeline.
//!
//! Mirrors STAR's intent (record-payload + iRead in the same buffer) but
//! uses a Rust-native layout instead of `BAMoutput.cpp`'s raw
//! `bam_bytes || uint iRead` concatenation. Phase 2b will port the native
//! layout alongside T1's byte-exact encoder.

use std::io::{Read, Write};

use anyhow::{bail, Context, Result};

/// Sentinel value written into `ref_id`/`pos` for unmapped records
/// (matches STAR's `bamIn32[1] == (uint32)-1` check in
/// `BAMoutput::coordOneAlign`).
pub const UNMAPPED_SENTINEL: u32 = u32::MAX;

/// Magic bytes identifying the phase-2a temp-record format ("RBR1").
pub const TMP_RECORD_MAGIC: u32 = 0x5242_5231;

/// Minimum record header size on disk (magic + len + ref_id + pos + i_read).
pub const TMP_RECORD_HEADER_BYTES: usize = 24;

/// One record as it lives inside a per-thread per-bin temp file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TmpRecord {
    pub ref_id: u32,
    pub pos: u32,
    pub i_read: u64,
    pub bam_bytes: Vec<u8>,
}

/// Random-access sort key for a record already materialised inside a bin.
///
/// `file_offset` is the byte offset (from the start of the loaded bin
/// buffer) where this record's header begins, so we can rewrite records
/// in sorted order without re-parsing each one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BinRecordMeta {
    pub coord: u64,
    pub i_read: u64,
    pub file_offset: u64,
}

impl BinRecordMeta {
    /// `(coord, i_read, file_offset)` — the STAR sort key, see
    /// `BAMbinSortByCoordinate.cpp:52` (`funCompareArrays<uint,3>`).
    pub fn sort_key(&self) -> (u64, u64, u64) {
        (self.coord, self.i_read, self.file_offset)
    }
}

/// Write one `TmpRecord` into a buffered writer. Little-endian on all
/// platforms so the format is portable across dev machines.
pub fn write_tmp_record<W: Write>(w: &mut W, rec: &TmpRecord) -> Result<()> {
    let rec_len: u32 = rec
        .bam_bytes
        .len()
        .try_into()
        .context("tmp record bam_bytes exceeds u32::MAX")?;
    w.write_all(&TMP_RECORD_MAGIC.to_le_bytes())?;
    w.write_all(&rec_len.to_le_bytes())?;
    w.write_all(&rec.ref_id.to_le_bytes())?;
    w.write_all(&rec.pos.to_le_bytes())?;
    w.write_all(&rec.i_read.to_le_bytes())?;
    w.write_all(&rec.bam_bytes)?;
    Ok(())
}

/// Read one `TmpRecord`. Returns `Ok(None)` at clean EOF, `Err` on
/// truncation or magic mismatch.
pub fn read_tmp_record<R: Read>(r: &mut R) -> Result<Option<TmpRecord>> {
    let mut header = [0u8; TMP_RECORD_HEADER_BYTES];
    let mut filled = 0usize;
    while filled < TMP_RECORD_HEADER_BYTES {
        let n = r.read(&mut header[filled..])?;
        if n == 0 {
            if filled == 0 {
                return Ok(None); // clean EOF
            }
            bail!("truncated tmp record header ({filled}/24 bytes)");
        }
        filled += n;
    }
    let magic = u32::from_le_bytes(header[0..4].try_into().unwrap());
    if magic != TMP_RECORD_MAGIC {
        bail!("tmp record magic mismatch: 0x{magic:08x}");
    }
    let rec_len = u32::from_le_bytes(header[4..8].try_into().unwrap()) as usize;
    let ref_id = u32::from_le_bytes(header[8..12].try_into().unwrap());
    let pos = u32::from_le_bytes(header[12..16].try_into().unwrap());
    let i_read = u64::from_le_bytes(header[16..24].try_into().unwrap());
    let mut bam_bytes = vec![0u8; rec_len];
    r.read_exact(&mut bam_bytes).context("reading tmp record payload")?;
    Ok(Some(TmpRecord { ref_id, pos, i_read, bam_bytes }))
}

/// Helper: 64-bit coordinate = (ref_id << 32) | pos, matching STAR's
/// `alignG = (bamIn32[1]<<32) | bamIn32[2]` packing in
/// `BAMoutput::coordOneAlign`. Returns `u64::MAX` for unmapped inputs.
pub fn pack_coord(ref_id: u32, pos: u32) -> u64 {
    if ref_id == UNMAPPED_SENTINEL {
        u64::MAX
    } else {
        ((ref_id as u64) << 32) | (pos as u64)
    }
}
```

- [ ] **Step 4: Register module**

Edit `crates/star-bam/src/lib.rs` — replace the existing module list with:

```rust
pub mod bam_cat;
pub mod bam_functions;
pub mod bam_output;
pub mod bam_remove_dups;
pub mod bam_sort_coord;
pub mod bam_sort_record;
pub mod sam_headers;
pub mod sam_output;
pub mod sam_to_bam;
pub mod signal_from_bam;
```

- [ ] **Step 5: Run tests to confirm pass**

Run: `cargo test -p star-bam --test bam_output_algorithm 2>&1 | tail -10`
Expected: `test result: ok. 4 passed` (4 tests from Step 1).

- [ ] **Step 6: Commit**

```bash
git add crates/star-bam/src/bam_sort_record.rs crates/star-bam/src/lib.rs crates/star-bam/tests/bam_output_algorithm.rs
git commit -m "star-bam: add BinRecordMeta + temp-record codec for phase-2a sort"
```

---

## Task 2: Parameters scaffolding

**Files:**
- Modify: `crates/star-params/src/parameters.rs:315-317` (add new fields), `:728` (defaults), around `:1247` (add CLI parsing).

- [ ] **Step 1: Write failing test**

Add to `crates/star-params/tests/parse.rs` (create if missing):

```rust
#[test]
fn parses_bam_sort_flags() {
    let mut p = star_params::parameters::Parameters::default_values();
    star_params::parameters::apply_args(
        &mut p,
        &[
            "--outBAMsortingThreadN", "6",
            "--outBAMsortingBinsN", "40",
            "--limitBAMsortRAM", "2000000000",
        ],
    ).unwrap();
    assert_eq!(p.out_bam_sorting_thread_n, 6);
    assert_eq!(p.out_bam_coord_nbins, 40);
    assert_eq!(p.limit_bam_sort_ram, 2_000_000_000);
}
```

Note: adjust to the existing test harness; if `apply_args` is named differently, mirror the test pattern already used in that file (find via `grep -rn "outBAMcompression" crates/star-params/tests/`).

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p star-params parses_bam_sort_flags 2>&1 | tail -20`
Expected: compile error — fields don't exist.

- [ ] **Step 3: Add fields**

In `crates/star-params/src/parameters.rs`, locate the struct field list after `pub out_bam_compression: i32,` (line 317) and insert:

```rust
    /// `--outBAMsortingBinsN` (STAR `outBAMcoordNbins`). Default 50.
    pub out_bam_coord_nbins: u32,
    /// `--outBAMsortingThreadN`. 0 = auto (use `run_thread_n`).
    pub out_bam_sorting_thread_n: u32,
    /// `--limitBAMsortRAM` in bytes. 0 = "use all available".
    pub limit_bam_sort_ram: u64,
    /// Runtime-populated genomic bin boundaries (length = `out_bam_coord_nbins`).
    /// Packed `(ref_id << 32) | pos`. Element 0 is always 0 once coord_bins
    /// runs. Populated by `BamOutput::coord_bins`.
    pub out_bam_sorting_bin_start: Vec<u64>,
    /// Directory for per-bin temp files during phase-2a sort (derived
    /// from `--outTmpDir`). Default: `"<outFileNamePrefix>_STARtmp/BAMsort"`.
    pub out_bam_sort_tmp_dir: String,
    /// Per-thread BAM chunk buffer budget (bytes). Default 16 MiB.
    pub chunk_out_bam_size_bytes: u64,
```

In the `default_values` impl (around line 728), add:

```rust
    out_bam_coord_nbins: 50,
    out_bam_sorting_thread_n: 0,
    limit_bam_sort_ram: 0,
    out_bam_sorting_bin_start: Vec::new(),
    out_bam_sort_tmp_dir: String::new(),
    chunk_out_bam_size_bytes: 16 * 1024 * 1024,
```

- [ ] **Step 4: Parse the flags**

Find the CLI arg match at `"outBAMcompression" => ...` (around line 1247) and add, in the same match block:

```rust
                    "outBAMsortingThreadN" => {
                        self.out_bam_sorting_thread_n = single()?.parse()?;
                    }
                    "outBAMsortingBinsN" => {
                        self.out_bam_coord_nbins = single()?.parse()?;
                    }
                    "limitBAMsortRAM" => {
                        self.limit_bam_sort_ram = single()?.parse()?;
                    }
```

- [ ] **Step 5: Run test**

Run: `cargo test -p star-params parses_bam_sort_flags 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add crates/star-params/src/parameters.rs crates/star-params/tests/parse.rs
git commit -m "star-params: add BAM-sort parameters (bins, threadN, limitRAM)"
```

---

## Task 3: `BamOutput::coord_one_align` (initial single-bin state)

**Files:**
- Modify: `crates/star-bam/src/bam_output.rs` (replace stub).
- Modify: `crates/star-bam/tests/bam_output_algorithm.rs` (add tests).
- Modify: `crates/star-bam/Cargo.toml` (add `noodles-bam` already present — confirm).

- [ ] **Step 1: Write failing tests**

Append to `crates/star-bam/tests/bam_output_algorithm.rs`:

```rust
use star_bam::bam_output::BamOutput;
use std::path::PathBuf;

fn tmp_dir(test_name: &str) -> PathBuf {
    let base = std::env::temp_dir().join(format!("star-bam-test-{test_name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&base);
    std::fs::create_dir_all(&base).unwrap();
    base
}

/// Build a fake BAM record bytes payload with the noodles record's
/// serialised form. For phase-2a tests we only care about ref_id/pos
/// extraction and the bam_bytes blob round-trip.
fn fake_bam_bytes(n: usize) -> Vec<u8> {
    // 4-byte block_size prefix (value unused by our sort path) + payload.
    let mut v = vec![0u8; 4 + n];
    v[0..4].copy_from_slice(&(n as u32).to_le_bytes());
    for (i, byte) in v[4..].iter_mut().enumerate() {
        *byte = (i % 251) as u8;
    }
    v
}

#[test]
fn coord_one_align_single_bin_accumulates() {
    let dir = tmp_dir("single-bin");
    let mut out = BamOutput::new(0, &dir, /* n_bins */ 50, /* chunk_bytes */ 1 << 16).unwrap();
    // feed 3 mapped records
    for i in 0..3u64 {
        out.coord_one_align(1, 1000 + i as u32, i, fake_bam_bytes(40)).unwrap();
    }
    // all should still be in bin 0 (pre-coord_bins)
    assert_eq!(out.bin_total_n()[0], 3);
    assert_eq!(out.active_bins(), 1);
}

#[test]
fn coord_one_align_unmapped_goes_to_last_bin() {
    let dir = tmp_dir("unmapped-last-bin");
    let mut out = BamOutput::new(0, &dir, /* n_bins */ 10, /* chunk_bytes */ 1 << 16).unwrap();
    out.coord_one_align(star_bam::bam_sort_record::UNMAPPED_SENTINEL, star_bam::bam_sort_record::UNMAPPED_SENTINEL, 99, fake_bam_bytes(40)).unwrap();
    // unmapped always lands in last bin (n_bins - 1 = 9), regardless of coord_bins state
    assert_eq!(out.bin_total_n()[9], 1);
    assert_eq!(out.bin_total_n()[0], 0);
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p star-bam --test bam_output_algorithm coord_one_align 2>&1 | head -30`
Expected: compile error — `BamOutput` struct not yet defined.

- [ ] **Step 3: Implement `BamOutput` with `coord_one_align` (coord_bins still a stub)**

Replace `crates/star-bam/src/bam_output.rs` entirely with:

```rust
//! Per-thread BAM record accumulator, mirroring `BAMoutput.cpp`'s coord
//! path. We keep STAR's semantics:
//!
//! * Start with `active_bins = 1`; everything except unmapped goes into
//!   bin 0 until `coord_bins` is called to pick real boundaries.
//! * Unmapped records always land in the last bin
//!   (`n_bins - 1`), matching STAR's `iBin = P.outBAMcoordNbins-1` path.
//! * Each bin has an in-memory buffer (flushed to a per-bin tmp file
//!   when full) plus totals (`bin_total_n`, `bin_total_bytes`) that the
//!   orchestrator reads to build the final assembly plan.
//!
//! Temp file layout per bin: see `bam_sort_record::TmpRecord`.

use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::bam_sort_record::{pack_coord, write_tmp_record, TmpRecord, UNMAPPED_SENTINEL};

pub struct BamOutput {
    /// Thread id / chunk index (`iChunk` in STAR). Phase 2a uses 0.
    pub i_chunk: u32,
    /// Directory containing this thread's per-bin temp files
    /// (`{tmp_root}/{i_chunk}/`).
    pub bin_dir: PathBuf,
    /// Configured number of bins (final, post-coord_bins value).
    pub n_bins_total: u32,
    /// Number of bins currently active. Starts at 1, jumps to
    /// `n_bins_total` after the first `coord_bins` call.
    pub n_bins_active: u32,
    /// Per-bin in-memory buffers (only used for mapped bins 0..n-2; bin
    /// n-1 and unmapped records currently bypass in-memory buffering
    /// for simplicity — they stream directly to their file writer).
    bin_buffers: Vec<Vec<u8>>,
    /// Per-bin file writers (one file per bin in this thread's dir).
    bin_writers: Vec<BufWriter<File>>,
    /// Per-bin capacity (bytes) for the in-memory buffer.
    bin_capacity: usize,
    /// STAR's `binSize1`: capacity while still in the single-bin state.
    bin_capacity_single: usize,
    /// Per-bin record count (STAR: `binTotalN`). Index 0..n_bins_total.
    bin_total_n: Vec<u64>,
    /// Per-bin total record-bytes written (STAR: `binTotalBytes`).
    bin_total_bytes: Vec<u64>,
    /// Cached records from the single-bin phase, needed by `coord_bins`
    /// to redistribute once real boundaries are known. Only populated
    /// while `n_bins_active == 1`.
    single_bin_cache: Vec<TmpRecord>,
    /// Runtime-computed bin boundaries (length n_bins_total - 1;
    /// element 0 is always 0). Populated by `coord_bins`.
    bin_starts: Vec<u64>,
}

impl BamOutput {
    pub fn new(i_chunk: u32, tmp_root: &Path, n_bins: u32, chunk_bytes: u64) -> Result<Self> {
        assert!(n_bins >= 2, "n_bins must be >= 2 (mapped bins + unmapped bin)");
        let bin_dir = tmp_root.join(i_chunk.to_string());
        create_dir_all(&bin_dir).with_context(|| format!("mkdir {}", bin_dir.display()))?;
        let per_bin_cap = (chunk_bytes / n_bins as u64) as usize;
        // STAR: `binSize1 = binStart[nBins-1] - binStart[0]`, i.e. the
        // entire mapped-bin allocation block minus the unmapped slot.
        let single_cap = per_bin_cap.saturating_mul((n_bins - 1) as usize);

        let mut bin_buffers = Vec::with_capacity(n_bins as usize);
        let mut bin_writers = Vec::with_capacity(n_bins as usize);
        for ib in 0..n_bins {
            bin_buffers.push(Vec::with_capacity(per_bin_cap));
            let fp = bin_dir.join(ib.to_string());
            let f = File::create(&fp).with_context(|| format!("create {}", fp.display()))?;
            bin_writers.push(BufWriter::with_capacity(1 << 16, f));
        }

        Ok(Self {
            i_chunk,
            bin_dir,
            n_bins_total: n_bins,
            n_bins_active: 1,
            bin_buffers,
            bin_writers,
            bin_capacity: per_bin_cap,
            bin_capacity_single: single_cap,
            bin_total_n: vec![0; n_bins as usize],
            bin_total_bytes: vec![0; n_bins as usize],
            single_bin_cache: Vec::new(),
            bin_starts: Vec::new(),
        })
    }

    /// Per-thread record count per bin (length = n_bins_total).
    pub fn bin_total_n(&self) -> &[u64] { &self.bin_total_n }
    /// Per-thread record-bytes per bin (length = n_bins_total).
    pub fn bin_total_bytes(&self) -> &[u64] { &self.bin_total_bytes }
    /// Number of active bins (1 until `coord_bins` runs).
    pub fn active_bins(&self) -> u32 { self.n_bins_active }

    /// Port of `BAMoutput::coordOneAlign`.
    ///
    /// `ref_id` and `pos` are extracted by the caller from the noodles
    /// record (use `UNMAPPED_SENTINEL` for both when the record is
    /// unmapped).
    pub fn coord_one_align(
        &mut self,
        ref_id: u32,
        pos: u32,
        i_read: u64,
        bam_bytes: Vec<u8>,
    ) -> Result<()> {
        if bam_bytes.is_empty() {
            return Ok(());
        }
        let is_unmapped = ref_id == UNMAPPED_SENTINEL;
        let i_bin = if is_unmapped {
            self.n_bins_total - 1
        } else if self.n_bins_active > 1 {
            self.pick_bin(pack_coord(ref_id, pos))
        } else {
            0
        };

        let rec = TmpRecord { ref_id, pos, i_read, bam_bytes };

        // Preserve single-bin-phase inputs verbatim for later redistribution.
        if self.n_bins_active == 1 && !is_unmapped {
            self.single_bin_cache.push(rec.clone());
        }

        self.push_to_bin(i_bin as usize, rec)
    }

    /// STAR's `binarySearch1a<uint64>(alignG, outBAMsortingBinStart, nBins-1)`
    /// picks the largest index `ib` such that `bin_starts[ib] <= alignG`.
    fn pick_bin(&self, align_g: u64) -> u32 {
        debug_assert!(!self.bin_starts.is_empty());
        let idx = match self.bin_starts.binary_search(&align_g) {
            Ok(i) => i,
            Err(i) => i.saturating_sub(1),
        };
        idx as u32
    }

    fn push_to_bin(&mut self, i_bin: usize, rec: TmpRecord) -> Result<()> {
        let payload_len = 24 + rec.bam_bytes.len();
        let cap = if self.n_bins_active > 1 || i_bin == (self.n_bins_total as usize - 1) {
            self.bin_capacity
        } else {
            self.bin_capacity_single
        };

        // Flush buffer if this record would overflow it.
        if self.bin_buffers[i_bin].len() + payload_len > cap {
            if self.n_bins_active > 1 || i_bin == (self.n_bins_total as usize - 1) {
                self.flush_bin_buffer(i_bin)?;
            } else {
                // STAR path: first-chunk overflow → run coord_bins and
                // redistribute, then re-attempt the current record.
                self.coord_bins()?;
                return self.coord_one_align(rec.ref_id, rec.pos, rec.i_read, rec.bam_bytes);
            }
        }

        // Append to in-memory bin buffer.
        let buf = &mut self.bin_buffers[i_bin];
        write_tmp_record(buf, &rec)?;
        self.bin_total_bytes[i_bin] += rec.bam_bytes.len() as u64;
        self.bin_total_n[i_bin] += 1;
        Ok(())
    }

    fn flush_bin_buffer(&mut self, i_bin: usize) -> Result<()> {
        let buf = &mut self.bin_buffers[i_bin];
        if buf.is_empty() {
            return Ok(());
        }
        self.bin_writers[i_bin].write_all(buf)?;
        buf.clear();
        Ok(())
    }

    /// Declared here so compile succeeds; fully implemented in Task 4.
    pub fn coord_bins(&mut self) -> Result<()> {
        // TASK 4 replaces this body.
        Ok(())
    }

    /// Declared here so compile succeeds; fully implemented in Task 5.
    pub fn coord_flush(&mut self) -> Result<()> {
        for ib in 0..self.n_bins_total as usize {
            self.flush_bin_buffer(ib)?;
            self.bin_writers[ib].flush()?;
        }
        Ok(())
    }
}
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p star-bam --test bam_output_algorithm coord_one_align 2>&1 | tail -15`
Expected: both `coord_one_align_single_bin_accumulates` and `coord_one_align_unmapped_goes_to_last_bin` PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/star-bam/src/bam_output.rs crates/star-bam/tests/bam_output_algorithm.rs
git commit -m "star-bam: implement BamOutput::coord_one_align (single-bin + unmapped-to-last)"
```

---

## Task 4: `BamOutput::coord_bins` boundary estimation

**Files:**
- Modify: `crates/star-bam/src/bam_output.rs` (replace `coord_bins` body).
- Modify: `crates/star-bam/tests/bam_output_algorithm.rs` (add tests).

- [ ] **Step 1: Write failing tests**

Append to `crates/star-bam/tests/bam_output_algorithm.rs`:

```rust
#[test]
fn coord_bins_splits_evenly_by_rank() {
    let dir = tmp_dir("coord-bins-even");
    // n_bins = 5 → 4 mapped bins + 1 unmapped bin.
    let mut out = BamOutput::new(0, &dir, 5, 1 << 16).unwrap();
    // Feed 8 mapped records across ref_id=1, positions 1000,2000,3000...8000.
    for i in 0..8u64 {
        out.coord_one_align(1, 1000 + 1000 * i as u32, i, fake_bam_bytes(32)).unwrap();
    }
    assert_eq!(out.active_bins(), 1);
    out.coord_bins().unwrap();
    assert_eq!(out.active_bins(), 5);
    // After coord_bins, bin 0 is the lowest-coord quarter (positions 1000,2000),
    // bin 1 is the next quarter, etc. Each of bins 0..=3 should now hold 2 records.
    let counts = out.bin_total_n();
    assert_eq!(counts[0..4].iter().sum::<u64>(), 8);
    // bin 4 is unmapped; no unmapped inputs → 0.
    assert_eq!(counts[4], 0);
}

#[test]
fn coord_bins_idempotent_after_first_call() {
    let dir = tmp_dir("coord-bins-idempotent");
    let mut out = BamOutput::new(0, &dir, 4, 1 << 16).unwrap();
    for i in 0..4u64 {
        out.coord_one_align(1, 100 * i as u32 + 100, i, fake_bam_bytes(24)).unwrap();
    }
    out.coord_bins().unwrap();
    let snapshot: Vec<u64> = out.bin_total_n().to_vec();
    out.coord_bins().unwrap();
    assert_eq!(out.bin_total_n(), snapshot.as_slice());
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p star-bam --test bam_output_algorithm coord_bins 2>&1 | tail -20`
Expected: tests compile but fail (`active_bins != 5`, `counts[0..4]` may be `[8,0,0,0]`).

- [ ] **Step 3: Implement `coord_bins`**

In `crates/star-bam/src/bam_output.rs`, replace the stub `coord_bins` with:

```rust
    /// Port of `BAMoutput::coordBins`.
    ///
    /// Called once, after enough records have been observed in bin 0 to
    /// estimate genomic boundaries. Reads the single-bin cache, extracts
    /// packed coords, sorts them, picks evenly-spaced rank cuts, stores
    /// them in `self.bin_starts`, then replays every cached record
    /// through `coord_one_align` so they land in their real bins.
    ///
    /// STAR reference: `BAMoutput.cpp:118-166`.
    pub fn coord_bins(&mut self) -> Result<()> {
        if self.n_bins_active != 1 {
            return Ok(()); // already done (STAR's mutex-guarded re-entry guard)
        }
        self.n_bins_active = self.n_bins_total;

        // Extract coords from cache and compute boundaries.
        let mut coords: Vec<u64> = self
            .single_bin_cache
            .iter()
            .map(|r| pack_coord(r.ref_id, r.pos))
            .collect();
        coords.sort_unstable();

        // STAR: `outBAMsortingBinStart[0] = 0;
        //        for ib in 1..(nBins-1):
        //            start[ib] = coords[ binTotalN[0] / (nBins-1) * ib ];`
        let mapped_bins = (self.n_bins_total - 1) as usize;
        let n = coords.len();
        let mut starts = vec![0u64; mapped_bins];
        for ib in 1..mapped_bins {
            let rank = n.saturating_mul(ib) / mapped_bins;
            let rank = rank.min(n.saturating_sub(1));
            starts[ib] = if n == 0 { 0 } else { coords[rank] };
        }
        self.bin_starts = starts;

        // Drain cache + per-bin totals for bin 0 (we're about to rewrite
        // them) and truncate bin 0's on-disk file: records that were
        // already spilled during the single-bin phase would double-count
        // otherwise. The in-memory buffer for bin 0 was where cached
        // records lived, so reset it too.
        let cached = std::mem::take(&mut self.single_bin_cache);
        self.bin_buffers[0].clear();
        self.bin_total_n[0] = 0;
        self.bin_total_bytes[0] = 0;
        // truncate bin-0 file to 0 bytes:
        self.bin_writers[0].flush()?;
        let path = self.bin_dir.join("0");
        {
            let f = File::create(&path)
                .with_context(|| format!("truncate {}", path.display()))?;
            self.bin_writers[0] = BufWriter::with_capacity(1 << 16, f);
        }

        // Redistribute every cached record.
        for rec in cached {
            self.coord_one_align(rec.ref_id, rec.pos, rec.i_read, rec.bam_bytes)?;
        }
        Ok(())
    }
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p star-bam --test bam_output_algorithm coord_bins 2>&1 | tail -15`
Expected: both `coord_bins_splits_evenly_by_rank` and `coord_bins_idempotent_after_first_call` PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/star-bam/src/bam_output.rs crates/star-bam/tests/bam_output_algorithm.rs
git commit -m "star-bam: implement BamOutput::coord_bins boundary estimation"
```

---

## Task 5: `BamOutput::coord_flush` + single-bin auto-coord_bins

**Files:**
- Modify: `crates/star-bam/src/bam_output.rs` (tighten `coord_flush`).
- Modify: `crates/star-bam/tests/bam_output_algorithm.rs` (add tests).

- [ ] **Step 1: Write failing test**

Append to `crates/star-bam/tests/bam_output_algorithm.rs`:

```rust
#[test]
fn coord_flush_without_explicit_coord_bins_still_distributes() {
    // STAR's coordFlush auto-calls coordBins if we're still in the
    // single-bin state — verify we do the same.
    let dir = tmp_dir("coord-flush-auto");
    let mut out = BamOutput::new(0, &dir, 4, 1 << 16).unwrap();
    for i in 0..6u64 {
        out.coord_one_align(1, 100 * (i as u32) + 10, i, fake_bam_bytes(32)).unwrap();
    }
    assert_eq!(out.active_bins(), 1);
    out.coord_flush().unwrap();
    assert_eq!(out.active_bins(), 4);
    let counts = out.bin_total_n();
    // mapped bins sum to 6; unmapped bin empty.
    assert_eq!(counts[0..3].iter().sum::<u64>(), 6);
    assert_eq!(counts[3], 0);
    // Per-thread bin files should exist and bin 0 should be non-empty.
    let bin0 = dir.join("0").join("0");
    let md = std::fs::metadata(&bin0).unwrap();
    assert!(md.len() > 0, "bin 0 should have been written to disk");
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p star-bam --test bam_output_algorithm coord_flush 2>&1 | tail -15`
Expected: `active_bins == 1` (auto-coord_bins not firing).

- [ ] **Step 3: Tighten `coord_flush`**

Replace `coord_flush` in `crates/star-bam/src/bam_output.rs`:

```rust
    /// Port of `BAMoutput::coordFlush`. If we never left the single-bin
    /// state, run `coord_bins` first so records are distributed before
    /// the final flush.
    pub fn coord_flush(&mut self) -> Result<()> {
        if self.n_bins_active == 1 {
            self.coord_bins()?;
        }
        for ib in 0..self.n_bins_total as usize {
            self.flush_bin_buffer(ib)?;
            self.bin_writers[ib].flush()?;
        }
        Ok(())
    }
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p star-bam --test bam_output_algorithm 2>&1 | tail -10`
Expected: all 7 tests PASS (4 from Task 1, 2 from Task 3, 2 from Task 4, 1 from this task — adjust if split differently).

- [ ] **Step 5: Commit**

```bash
git add crates/star-bam/src/bam_output.rs crates/star-bam/tests/bam_output_algorithm.rs
git commit -m "star-bam: coord_flush auto-invokes coord_bins when still single-binned"
```

---

## Task 6: `sort_mapped_bin` per-bin coordinate sort

**Files:**
- Create: `crates/star-bam/src/bam_sort_bin.rs`.
- Modify: `crates/star-bam/src/lib.rs` (register module).
- Modify: `crates/star-bam/tests/bam_output_algorithm.rs` (add test).

- [ ] **Step 1: Write failing test**

Append to `crates/star-bam/tests/bam_output_algorithm.rs`:

```rust
use star_bam::bam_sort_bin::sort_mapped_bin;
use star_bam::bam_sort_record::read_tmp_record;
use std::fs::File;
use std::io::BufReader;

fn write_tmp_bin_file(path: &std::path::Path, recs: &[star_bam::bam_sort_record::TmpRecord]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for r in recs {
        star_bam::bam_sort_record::write_tmp_record(&mut f, r).unwrap();
    }
}

#[test]
fn sort_mapped_bin_sorts_by_coord_then_iread() {
    let root = tmp_dir("sort-mapped-bin");
    // Two threads, one bin (ibin=0).
    for t in 0..2u32 {
        std::fs::create_dir_all(root.join(t.to_string())).unwrap();
    }
    use star_bam::bam_sort_record::TmpRecord;
    let bam = |n: usize| fake_bam_bytes(n);
    write_tmp_bin_file(&root.join("0").join("0"), &[
        TmpRecord { ref_id: 1, pos: 300, i_read: 3, bam_bytes: bam(20) },
        TmpRecord { ref_id: 1, pos: 100, i_read: 1, bam_bytes: bam(20) },
    ]);
    write_tmp_bin_file(&root.join("1").join("0"), &[
        TmpRecord { ref_id: 1, pos: 200, i_read: 2, bam_bytes: bam(20) },
        TmpRecord { ref_id: 1, pos: 100, i_read: 0, bam_bytes: bam(20) },
    ]);

    let bin_n = 4u64;
    let bin_s = 4 * (24 + 24); // 4 records × (24 header + 24 payload incl prefix)
    sort_mapped_bin(0, bin_n, bin_s as u64, 2, &root).unwrap();

    // Read back the sorted output file and check i_read order.
    let out_path = root.join("bin_0000_sorted.rbr");
    let f = File::open(&out_path).unwrap();
    let mut r = BufReader::new(f);
    let mut iread_order = Vec::new();
    while let Some(rec) = read_tmp_record(&mut r).unwrap() {
        iread_order.push(rec.i_read);
    }
    // Coord 100 ties (iread 0, 1), then 200 (iread 2), then 300 (iread 3).
    assert_eq!(iread_order, vec![0, 1, 2, 3]);
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p star-bam --test bam_output_algorithm sort_mapped_bin 2>&1 | tail -20`
Expected: compile error — `bam_sort_bin` module missing.

- [ ] **Step 3: Create `bam_sort_bin.rs`**

Create `crates/star-bam/src/bam_sort_bin.rs`:

```rust
//! Per-bin mapped-record sort — phase-2a port of
//! `BAMbinSortByCoordinate.cpp`.
//!
//! For one bin index, load every per-thread tmp file for that bin into a
//! single flat buffer, build a `BinRecordMeta` array keyed by
//! `(coord, i_read, offset)`, sort, and write the records back out in
//! sorted order. Output file: `{tmp_root}/bin_{ibin:04}_sorted.rbr`.
//!
//! Output format is the same `TmpRecord` codec used everywhere in the
//! phase-2a pipeline — the final-assembly step (Task 9) decodes it and
//! writes the BAM via `noodles-bam`.

use std::fs::{remove_file, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::bam_sort_record::{
    pack_coord, read_tmp_record, write_tmp_record, BinRecordMeta, TmpRecord,
};

/// Port of `BAMbinSortByCoordinate`.
///
/// * `ibin` — the bin index (0..n_bins_total-2 for mapped bins).
/// * `bin_n` — expected record count (for a sanity check only).
/// * `bin_s` — expected total bytes across all thread files (sanity check).
/// * `n_threads` — number of per-thread input files to read.
/// * `tmp_root` — the parent dir containing per-thread subdirs.
pub fn sort_mapped_bin(
    ibin: u32,
    bin_n: u64,
    bin_s: u64,
    n_threads: u32,
    tmp_root: &Path,
) -> Result<()> {
    if bin_n == 0 {
        return Ok(());
    }

    // Load every per-thread tmp file for this bin into one flat buffer.
    let mut all: Vec<u8> = Vec::with_capacity(bin_s as usize);
    for it in 0..n_threads {
        let path = tmp_root.join(it.to_string()).join(ibin.to_string());
        if !path.exists() {
            continue;
        }
        let mut f = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;
        all.extend_from_slice(&buf);
    }

    // Parse records into `(BinRecordMeta, TmpRecord)`. We keep the
    // TmpRecord around so we can emit in sorted order without a second
    // parse pass; file_offset in `BinRecordMeta` becomes the index into
    // the parallel Vec, not a byte offset (equivalent for sorting).
    let mut records: Vec<TmpRecord> = Vec::with_capacity(bin_n as usize);
    let mut metas: Vec<BinRecordMeta> = Vec::with_capacity(bin_n as usize);
    let mut cursor = std::io::Cursor::new(&all[..]);
    let mut idx: u64 = 0;
    while let Some(rec) = read_tmp_record(&mut cursor)? {
        let coord = pack_coord(rec.ref_id, rec.pos);
        metas.push(BinRecordMeta { coord, i_read: rec.i_read, file_offset: idx });
        records.push(rec);
        idx += 1;
    }
    if metas.len() as u64 != bin_n {
        bail!(
            "bin {ibin}: expected {bin_n} records, found {} on disk",
            metas.len()
        );
    }

    metas.sort_by(|a, b| a.sort_key().cmp(&b.sort_key()));

    // Emit in sorted order.
    let out_path = tmp_root.join(format!("bin_{ibin:04}_sorted.rbr"));
    let out = File::create(&out_path).with_context(|| format!("create {}", out_path.display()))?;
    let mut w = BufWriter::with_capacity(1 << 16, out);
    for m in &metas {
        write_tmp_record(&mut w, &records[m.file_offset as usize])?;
    }
    w.flush()?;

    // Clean up per-thread source files.
    for it in 0..n_threads {
        let path = tmp_root.join(it.to_string()).join(ibin.to_string());
        if path.exists() {
            let _ = remove_file(&path);
        }
    }
    Ok(())
}
```

- [ ] **Step 4: Register module**

Edit `crates/star-bam/src/lib.rs` module list — insert after `bam_sort_coord`:

```rust
pub mod bam_sort_bin;
```

- [ ] **Step 5: Run test**

Run: `cargo test -p star-bam --test bam_output_algorithm sort_mapped_bin 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add crates/star-bam/src/bam_sort_bin.rs crates/star-bam/src/lib.rs crates/star-bam/tests/bam_output_algorithm.rs
git commit -m "star-bam: sort_mapped_bin — per-bin coord+iRead sort"
```

---

## Task 7: `merge_unmapped_bin` multi-way iRead merge

**Files:**
- Create: `crates/star-bam/src/bam_sort_unmapped.rs`.
- Modify: `crates/star-bam/src/lib.rs` (register module).
- Modify: `crates/star-bam/tests/bam_output_algorithm.rs` (add test).

- [ ] **Step 1: Write failing test**

Append to `crates/star-bam/tests/bam_output_algorithm.rs`:

```rust
use star_bam::bam_sort_unmapped::merge_unmapped_bin;

#[test]
fn merge_unmapped_bin_orders_by_iread_across_threads() {
    let root = tmp_dir("merge-unmapped");
    // Two threads, unmapped bin index = 3 (n_bins_total = 4).
    use star_bam::bam_sort_record::{TmpRecord, UNMAPPED_SENTINEL};
    for t in 0..2u32 {
        std::fs::create_dir_all(root.join(t.to_string())).unwrap();
    }
    // Per-thread inputs are already in ascending iRead order (what
    // coord_flush gives us because records are fed in input order).
    write_tmp_bin_file(&root.join("0").join("3"), &[
        TmpRecord { ref_id: UNMAPPED_SENTINEL, pos: UNMAPPED_SENTINEL, i_read: 2, bam_bytes: fake_bam_bytes(16) },
        TmpRecord { ref_id: UNMAPPED_SENTINEL, pos: UNMAPPED_SENTINEL, i_read: 5, bam_bytes: fake_bam_bytes(16) },
    ]);
    write_tmp_bin_file(&root.join("1").join("3"), &[
        TmpRecord { ref_id: UNMAPPED_SENTINEL, pos: UNMAPPED_SENTINEL, i_read: 1, bam_bytes: fake_bam_bytes(16) },
        TmpRecord { ref_id: UNMAPPED_SENTINEL, pos: UNMAPPED_SENTINEL, i_read: 3, bam_bytes: fake_bam_bytes(16) },
        TmpRecord { ref_id: UNMAPPED_SENTINEL, pos: UNMAPPED_SENTINEL, i_read: 4, bam_bytes: fake_bam_bytes(16) },
    ]);

    merge_unmapped_bin(3, 2, &root).unwrap();

    let out_path = root.join("bin_0003_sorted.rbr");
    let f = std::fs::File::open(&out_path).unwrap();
    let mut r = std::io::BufReader::new(f);
    let mut order = Vec::new();
    while let Some(rec) = star_bam::bam_sort_record::read_tmp_record(&mut r).unwrap() {
        order.push(rec.i_read);
    }
    assert_eq!(order, vec![1, 2, 3, 4, 5]);
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p star-bam --test bam_output_algorithm merge_unmapped_bin 2>&1 | tail -20`
Expected: compile error — `bam_sort_unmapped` missing.

- [ ] **Step 3: Create `bam_sort_unmapped.rs`**

Create `crates/star-bam/src/bam_sort_unmapped.rs`:

```rust
//! Phase-2a port of `BAMbinSortUnmapped.cpp`.
//!
//! Each per-thread unmapped bin temp file is already in ascending `iRead`
//! order (records were fed to `BamOutput::coord_one_align` in input
//! order within each thread). We perform a classic k-way merge using a
//! `BinaryHeap<Reverse<_>>` keyed on `i_read`.
//!
//! Output: `{tmp_root}/bin_{ibin:04}_sorted.rbr` — same `TmpRecord`
//! codec as the mapped bins, so the final assembly step can concat them
//! with zero format plumbing.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{remove_file, File};
use std::io::{BufReader, BufWriter};
use std::path::Path;

use anyhow::{Context, Result};

use crate::bam_sort_record::{read_tmp_record, write_tmp_record, TmpRecord};

struct HeapEntry {
    i_read: u64,
    stream_idx: usize,
    rec: TmpRecord,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool { self.i_read == other.i_read && self.stream_idx == other.stream_idx }
}
impl Eq for HeapEntry {}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> { Some(self.cmp(other)) }
}
impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Smaller i_read first, then smaller stream_idx to break ties
        // deterministically.
        self.i_read
            .cmp(&other.i_read)
            .then(self.stream_idx.cmp(&other.stream_idx))
    }
}

pub fn merge_unmapped_bin(ibin: u32, n_threads: u32, tmp_root: &Path) -> Result<()> {
    // Open each per-thread file that exists and has at least one record.
    let mut streams: Vec<BufReader<File>> = Vec::with_capacity(n_threads as usize);
    let mut paths: Vec<std::path::PathBuf> = Vec::new();
    let mut heap: BinaryHeap<Reverse<HeapEntry>> = BinaryHeap::new();

    for it in 0..n_threads {
        let path = tmp_root.join(it.to_string()).join(ibin.to_string());
        if !path.exists() {
            continue;
        }
        let f = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let mut r = BufReader::with_capacity(1 << 16, f);
        if let Some(first) = read_tmp_record(&mut r)? {
            let idx = streams.len();
            heap.push(Reverse(HeapEntry { i_read: first.i_read, stream_idx: idx, rec: first }));
            streams.push(r);
            paths.push(path);
        }
    }

    let out_path = tmp_root.join(format!("bin_{ibin:04}_sorted.rbr"));
    let out = File::create(&out_path).with_context(|| format!("create {}", out_path.display()))?;
    let mut w = BufWriter::with_capacity(1 << 16, out);

    while let Some(Reverse(entry)) = heap.pop() {
        write_tmp_record(&mut w, &entry.rec)?;
        // Pull the next record from the same stream.
        if let Some(next) = read_tmp_record(&mut streams[entry.stream_idx])? {
            heap.push(Reverse(HeapEntry { i_read: next.i_read, stream_idx: entry.stream_idx, rec: next }));
        }
    }
    drop(w);

    for p in &paths {
        let _ = remove_file(p);
    }
    Ok(())
}
```

- [ ] **Step 4: Register module**

Edit `crates/star-bam/src/lib.rs` — add:

```rust
pub mod bam_sort_unmapped;
```

- [ ] **Step 5: Run test**

Run: `cargo test -p star-bam --test bam_output_algorithm merge_unmapped_bin 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add crates/star-bam/src/bam_sort_unmapped.rs crates/star-bam/src/lib.rs crates/star-bam/tests/bam_output_algorithm.rs
git commit -m "star-bam: merge_unmapped_bin — k-way iRead merge"
```

---

## Task 8: `feed_sam_into_bam_output` (phase-2a SAM → BamOutput adapter)

**Files:**
- Create: `crates/star-bam/src/bam_sort_feed.rs`.
- Modify: `crates/star-bam/src/lib.rs` (register module).
- Modify: `crates/star-bam/tests/bam_output_algorithm.rs` (add test).

- [ ] **Step 1: Write failing test**

Append to `crates/star-bam/tests/bam_output_algorithm.rs`:

```rust
use star_bam::bam_sort_feed::feed_sam_into_bam_output;

#[test]
fn feed_sam_into_bam_output_assigns_monotonic_iread() {
    // Hand-build a tiny SAM with header + 3 records: one mapped to ref 0,
    // one mapped to ref 0 at a later pos, one unmapped.
    let sam = "\
@HD\tVN:1.6\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:1000\n\
r1\t0\tchr1\t10\t30\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n\
r2\t0\tchr1\t50\t30\t10M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII\n\
r3\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGGGGGGG\tIIIIIIIIII\n";
    let dir = tmp_dir("feed-sam");
    let sam_path = dir.join("in.sam");
    std::fs::write(&sam_path, sam).unwrap();

    let mut out = BamOutput::new(0, &dir.join("bamtmp"), 4, 1 << 16).unwrap();
    feed_sam_into_bam_output(&sam_path, &mut out).unwrap();
    out.coord_flush().unwrap();

    // 2 mapped (bins 0..2) + 1 unmapped (bin 3); iReads 0,1,2 assigned.
    let counts = out.bin_total_n();
    assert_eq!(counts.iter().sum::<u64>(), 3);
    assert_eq!(counts[3], 1, "unmapped lands in last bin");
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p star-bam --test bam_output_algorithm feed_sam_into 2>&1 | tail -20`
Expected: compile error — `bam_sort_feed` missing.

- [ ] **Step 3: Create `bam_sort_feed.rs`**

Create `crates/star-bam/src/bam_sort_feed.rs`:

```rust
//! Phase-2a adapter: read the existing `Aligned.out.sam.tmp` SAM
//! intermediate, convert every record to `noodles_bam` wire bytes, and
//! feed them to a `BamOutput` coord pipeline.
//!
//! `iRead` is assigned as a monotonic counter starting at 0 — the
//! "stable within thread" guarantee (phase 2a does not try to match
//! C++ STAR's global-ordinal iRead).

use std::fs::File;
use std::io::{BufReader, Cursor};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam as sam;

use crate::bam_output::BamOutput;
use crate::bam_sort_record::UNMAPPED_SENTINEL;

/// Drive a full SAM intermediate through `BamOutput::coord_one_align`.
///
/// Returns the number of records fed (for logging / stats).
pub fn feed_sam_into_bam_output(sam_path: &Path, out: &mut BamOutput) -> Result<u64> {
    let file = File::open(sam_path).with_context(|| format!("open {}", sam_path.display()))?;
    let mut reader = sam::io::Reader::new(BufReader::new(file));
    let header = reader.read_header().context("parse SAM header")?;

    let mut sam_rec = sam::Record::default();
    let mut i_read: u64 = 0;
    let mut n_fed: u64 = 0;
    loop {
        let n = reader.read_record(&mut sam_rec).context("reading SAM record")?;
        if n == 0 { break; }

        // Serialise this record as BAM wire bytes via a throwaway
        // in-memory BAM writer: this is the phase-2a hack. Phase 2b
        // swaps this for records produced natively by the ported
        // alignBAM encoder.
        let mut buf: Vec<u8> = Vec::with_capacity(256);
        {
            let mut bw = bam::io::Writer::from(Cursor::new(&mut buf));
            // Write a one-off record directly (no header — we want just
            // the body bytes).
            use sam::alignment::io::Write as _;
            bw.write_alignment_record(&header, &sam_rec)?;
        }
        // `noodles` prepends a BGZF + BAM header when instantiated via
        // Writer::from; strip header/BGZF framing so `buf` contains
        // just the record's serialised BAM bytes.
        //
        // Simpler path: call `bam::record::codec::encoder::encode_record`
        // directly — but it's not public. Fall back to encoding via
        // `bam::Record::try_from(&sam::Record)` and writing that.
        //
        // The ergonomic path: decode the SAM record into a
        // `bam::Record`, then write via `bam::io::Writer` into a
        // `Vec<u8>` ignoring the BAM header framing — but noodles'
        // Writer wraps in BGZF. Instead we use `bam::Record::try_from`.
        let bam_rec = bam::Record::try_from_alignment_record(&header, &sam_rec)
            .context("sam::Record → bam::Record")?;
        // `bam::Record` is already the on-disk BAM record (block_size
        // + body); serialise via its own `as_bytes`-equivalent. noodles
        // exposes this through `bam::io::Writer::write_record` on a
        // cursor; since that also doesn't emit header-less output, the
        // simplest reliable approach is to re-implement the wire
        // layout here via `encoder::record::encode_record` in noodles.
        // To avoid chasing private APIs, we sidestep: the phase-2a
        // pipeline never inspects `bam_bytes` field-by-field — it only
        // needs ref_id / pos (derived from `bam_rec.reference_sequence_id`
        // / `bam_rec.alignment_start`) and a roundtrippable blob.
        //
        // So: store the serialised `bam_rec` via `bincode` or use its
        // `as_ref()` bytes. noodles-bam 0.88's `bam::Record` stores
        // the on-disk representation; its `Deref<Target=[u8]>` returns
        // the record bytes with block_size prefix. Use that.
        let bam_bytes: Vec<u8> = bam_rec.as_ref().to_vec();

        let (ref_id, pos) = if let Some(rid) = bam_rec.reference_sequence_id().transpose()? {
            let rid = rid as u32;
            let pos = bam_rec
                .alignment_start()
                .transpose()?
                .map(|p| (usize::from(p) - 1) as u32)
                .unwrap_or(0);
            (rid, pos)
        } else {
            (UNMAPPED_SENTINEL, UNMAPPED_SENTINEL)
        };

        out.coord_one_align(ref_id, pos, i_read, bam_bytes)?;
        i_read += 1;
        n_fed += 1;
    }
    Ok(n_fed)
}
```

**⚠️ Implementation note for the agent executing this task.** The noodles API surface for "give me just the BAM-wire bytes of one record" has moved between 0.85 and 0.88. Before writing the code above, verify the two calls used:

1. `bam::Record::try_from_alignment_record(&header, &sam_rec)` — if not present under that exact name, search noodles-bam 0.88's docs for the `TryFrom<&sam::Record>` impl or the `bam::record::Record::from_alignment_record` free function. The semantics you need are "sam::Record + header → bam::Record in on-disk byte form."
2. `bam_rec.as_ref()` returning `&[u8]` — confirm by running `cargo doc -p noodles-bam --no-deps --open` locally and checking the `Record` docs. If `AsRef<[u8]>` isn't there, look for `.as_bytes()`, `.into_bytes()`, or the `bam::record::Record` field that holds the serialized form.

If the API diverges, adapt — but keep the contract identical: given a `sam::Record + Header`, produce `(ref_id, pos, bam_bytes)` where `bam_bytes` is the record's on-disk BAM serialisation including the 4-byte block_size prefix. Then feed that to `coord_one_align`. Do not change `coord_one_align`'s signature.

- [ ] **Step 4: Register module**

Edit `crates/star-bam/src/lib.rs`:

```rust
pub mod bam_sort_feed;
```

- [ ] **Step 5: Run tests**

Run: `cargo test -p star-bam --test bam_output_algorithm feed_sam 2>&1 | tail -15`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add crates/star-bam/src/bam_sort_feed.rs crates/star-bam/src/lib.rs crates/star-bam/tests/bam_output_algorithm.rs
git commit -m "star-bam: feed_sam_into_bam_output — phase-2a SAM→BAMoutput adapter"
```

---

## Task 9: `sort_bam_by_coordinate` orchestrator + final assembly

**Files:**
- Modify: `crates/star-bam/src/bam_sort_coord.rs` (replace stub).
- Modify: `crates/star-bam/Cargo.toml` (add `rayon.workspace = true` if not present).
- Modify: `crates/star-bam/tests/bam_output_algorithm.rs` (add test).

- [ ] **Step 1: Add rayon to star-bam**

Edit `crates/star-bam/Cargo.toml` — add under `[dependencies]`:

```toml
rayon = { workspace = true }
```

- [ ] **Step 2: Write failing test**

Append to `crates/star-bam/tests/bam_output_algorithm.rs`:

```rust
use star_bam::bam_sort_coord::sort_bam_by_coordinate;

#[test]
fn sort_bam_by_coordinate_end_to_end_tiny() {
    // Full pipeline: one thread, one BamOutput, feed via coord_one_align,
    // run sort, check final BAM exists and contains the right # of records.
    let dir = tmp_dir("sort-bam-e2e");
    let tmp_root = dir.join("bamsort");
    std::fs::create_dir_all(&tmp_root).unwrap();
    // Minimal SAM header for the output BAM.
    let sam_hdr = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";

    let mut out = BamOutput::new(0, &tmp_root, 4, 1 << 16).unwrap();
    // 3 mapped (coord ascending already) + 1 unmapped.
    let sam = "\
@HD\tVN:1.6\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:1000\n\
r1\t0\tchr1\t10\t30\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n\
r2\t0\tchr1\t50\t30\t10M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII\n\
r3\t0\tchr1\t30\t30\t10M\t*\t0\t0\tGGGGGGGGGG\tIIIIIIIIII\n\
r4\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\n";
    let sam_path = dir.join("in.sam");
    std::fs::write(&sam_path, sam).unwrap();
    feed_sam_into_bam_output(&sam_path, &mut out).unwrap();
    out.coord_flush().unwrap();

    let final_bam = dir.join("Aligned.sortedByCoord.out.bam");
    sort_bam_by_coordinate(
        &tmp_root,
        /* n_threads */ 1,
        /* n_bins_total */ 4,
        &[out.bin_total_n().to_vec()],
        &[out.bin_total_bytes().to_vec()],
        sam_hdr.as_bytes(),
        &final_bam,
        /* compression */ Some(1),
        /* limit_bam_sort_ram */ 0, // 0 = no limit in phase-2a
    )
    .unwrap();

    assert!(final_bam.exists());
    // Open final BAM and count records.
    let f = std::fs::File::open(&final_bam).unwrap();
    let mut r = noodles_bam::io::Reader::new(f);
    let _hdr = r.read_header().unwrap();
    let n: usize = r.records().count();
    assert_eq!(n, 4, "final BAM should contain all 4 records");
}
```

Add the import `use noodles_bam;` at the top if not yet in scope.

- [ ] **Step 3: Run to confirm failure**

Run: `cargo test -p star-bam --test bam_output_algorithm sort_bam_by_coordinate 2>&1 | tail -20`
Expected: compile error — `sort_bam_by_coordinate` not yet implemented with this signature.

- [ ] **Step 4: Implement orchestrator**

Replace `crates/star-bam/src/bam_sort_coord.rs` entirely with:

```rust
//! Phase-2a port of `bamSortByCoordinate.cpp`.
//!
//! Responsibilities:
//! 1. Aggregate per-thread `binTotalN` / `binTotalBytes` stats.
//! 2. Compute peak-memory estimate; bail if it exceeds
//!    `limit_bam_sort_ram` (skipped when limit == 0).
//! 3. Run per-bin sort in parallel via rayon:
//!    * mapped bins → `sort_mapped_bin`
//!    * unmapped bin (last index) → `merge_unmapped_bin`
//! 4. Assemble the final BAM by reading each per-bin `bin_{ibin:04}_sorted.rbr`
//!    in ascending bin order and streaming records through
//!    `noodles_bam::io::Writer`.

use std::fs::File;
use std::io::{BufReader, BufWriter, Cursor, Read};
use std::path::Path;

use anyhow::{bail, Context, Result};
use noodles_bam as bam;
use noodles_sam as sam;
use rayon::prelude::*;

use crate::bam_sort_bin::sort_mapped_bin;
use crate::bam_sort_record::read_tmp_record;
use crate::bam_sort_unmapped::merge_unmapped_bin;

/// Collect per-bin size / count across all threads.
pub fn collect_bin_stats(
    per_thread_n: &[Vec<u64>],
    per_thread_bytes: &[Vec<u64>],
    n_bins_total: u32,
) -> (Vec<u64>, Vec<u64>) {
    let mut n = vec![0u64; n_bins_total as usize];
    let mut b = vec![0u64; n_bins_total as usize];
    for t in 0..per_thread_n.len() {
        for ib in 0..n_bins_total as usize {
            n[ib] += per_thread_n[t][ib];
            b[ib] += per_thread_bytes[t][ib];
        }
    }
    (n, b)
}

/// STAR: `if (maxMem > limitBAMsortRAM) exit`. Returns the estimated
/// peak mapped-bin memory in bytes. Skips the check when `limit == 0`.
pub fn check_sort_ram_limit(
    per_bin_n: &[u64],
    per_bin_bytes: &[u64],
    limit: u64,
) -> Result<u64> {
    // STAR formula: per-bin mem = bytes + 24 * n (the 24 matches the
    // BAMoutput tmp-record header size; close enough for our Rust codec).
    let mapped = per_bin_n.len().saturating_sub(1);
    let mut peak = 0u64;
    for ib in 0..mapped {
        let m = per_bin_bytes[ib] + 24 * per_bin_n[ib];
        if m > peak { peak = m; }
    }
    if limit > 0 && peak > limit {
        bail!(
            "BAM sort would need {peak} bytes (> --limitBAMsortRAM={limit}). \
             Re-run with --limitBAMsortRAM {}",
            peak + 1_000_000_000
        );
    }
    Ok(peak)
}

#[allow(clippy::too_many_arguments)]
pub fn sort_bam_by_coordinate(
    tmp_root: &Path,
    n_threads: u32,
    n_bins_total: u32,
    per_thread_n: &[Vec<u64>],
    per_thread_bytes: &[Vec<u64>],
    sam_header_bytes: &[u8],
    final_bam_path: &Path,
    compression: Option<u32>,
    limit_bam_sort_ram: u64,
) -> Result<()> {
    let (bin_n, bin_s) = collect_bin_stats(per_thread_n, per_thread_bytes, n_bins_total);
    check_sort_ram_limit(&bin_n, &bin_s, limit_bam_sort_ram)?;

    let unmapped_ix = (n_bins_total - 1) as usize;
    // Parallel fan-out: mapped bins 0..unmapped_ix + the unmapped bin.
    // Process unmapped first (largest, and its scheduling-order in STAR
    // is "reverse" — unmapped bin first).
    let all_bins: Vec<u32> = (0..n_bins_total).rev().collect();
    all_bins.par_iter().try_for_each(|&ibin| -> Result<()> {
        let i = ibin as usize;
        if bin_n[i] == 0 && ibin != (n_bins_total - 1) {
            return Ok(());
        }
        if ibin == (n_bins_total - 1) {
            merge_unmapped_bin(ibin, n_threads, tmp_root)?;
        } else {
            sort_mapped_bin(ibin, bin_n[i], bin_s[i], n_threads, tmp_root)?;
        }
        Ok(())
    })?;

    assemble_final_bam(
        tmp_root,
        n_bins_total,
        sam_header_bytes,
        final_bam_path,
        compression,
    )
}

/// Read the SAM header supplied by the caller, decode it via
/// noodles-sam, then stream every per-bin `bin_{ibin:04}_sorted.rbr` in
/// ascending bin order, decoding each TmpRecord back into a
/// `bam::Record` and writing it through `bam::io::Writer`.
pub fn assemble_final_bam(
    tmp_root: &Path,
    n_bins_total: u32,
    sam_header_bytes: &[u8],
    final_bam_path: &Path,
    compression: Option<u32>,
) -> Result<()> {
    // Parse header.
    let mut hdr_reader = sam::io::Reader::new(Cursor::new(sam_header_bytes));
    let mut header = hdr_reader.read_header().context("parse SAM header for final BAM")?;
    // STAR writes "SO:coordinate" in samHeaderSortedCoord; we patch the
    // header here so the output SO matches.
    use sam::header::record::value::map::header::SortOrder;
    if let Some(h) = header.header_mut() {
        h.set_sort_order(SortOrder::Coordinate);
    }

    let out_file = File::create(final_bam_path)
        .with_context(|| format!("create {}", final_bam_path.display()))?;
    let bgzf_writer = if let Some(level) = compression {
        let lvl = noodles_bgzf::io::writer::CompressionLevel::try_from(level as u8)
            .unwrap_or_default();
        noodles_bgzf::io::writer::Builder::default()
            .set_compression_level(lvl)
            .build_from_writer(BufWriter::new(out_file))
    } else {
        noodles_bgzf::io::Writer::new(BufWriter::new(out_file))
    };
    let mut bam_w = bam::io::Writer::from(bgzf_writer);
    bam_w.write_header(&header)?;

    for ibin in 0..n_bins_total {
        let path = tmp_root.join(format!("bin_{ibin:04}_sorted.rbr"));
        if !path.exists() {
            continue;
        }
        let f = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let mut r = BufReader::with_capacity(1 << 16, f);
        while let Some(rec) = read_tmp_record(&mut r)? {
            // Reconstruct bam::Record from stored bytes.
            let bam_rec: bam::Record =
                bam::Record::try_from(rec.bam_bytes.as_slice()).with_context(|| {
                    format!("decoding TmpRecord → bam::Record (bin {ibin}, iRead {})", rec.i_read)
                })?;
            bam_w.write_record(&header, &bam_rec)?;
        }
        let _ = std::fs::remove_file(&path);
    }
    bam_w.try_finish()?;
    Ok(())
}
```

**⚠️ Same API-drift caveat as Task 8, step 3.** `bam::Record::try_from(&[u8])` and `bam_w.write_record(&header, &bam_rec)` are the phase-2a contract. If the noodles-bam 0.88 API names differ (e.g. `write_record` wants a different signature, or `try_from` needs a BGZF-framed input), adapt while keeping the contract: "given a `Vec<u8>` previously obtained from a `bam::Record`, round-trip it back into the writer."

- [ ] **Step 5: Run tests**

Run: `cargo test -p star-bam --test bam_output_algorithm sort_bam 2>&1 | tail -15`
Expected: `sort_bam_by_coordinate_end_to_end_tiny` PASS.

- [ ] **Step 6: Commit**

```bash
git add crates/star-bam/src/bam_sort_coord.rs crates/star-bam/Cargo.toml crates/star-bam/tests/bam_output_algorithm.rs
git commit -m "star-bam: sort_bam_by_coordinate orchestrator + final assembly"
```

---

## Task 10: Wire SortedByCoordinate path in `star-cli`

**Files:**
- Modify: `crates/star-cli/src/main.rs:259-263, 394-405`.

- [ ] **Step 1: Read the current branch logic**

Run: `rg -n "want_bam|out_bam_coord|sam_scratch_path" crates/star-cli/src/main.rs`
Confirm the branch is still at the lines above before editing.

- [ ] **Step 2: Update fallback comment + split branch**

Edit `crates/star-cli/src/main.rs` around lines 257-263. Replace:

```rust
    // Decide whether to keep the SAM on disk (`--outSAMtype SAM` or nothing
    // set — default), or to transiently produce SAM and then
    // convert to BAM (`--outSAMtype BAM Unsorted`). `SortedByCoordinate`
    // is not yet implemented — we currently fall back to Unsorted for
    // any BAM request.
    let want_sam = emit_sam && (p.out_sam_bool || (!p.out_bam_unsorted && !p.out_bam_coord));
    let want_bam = emit_sam && (p.out_bam_unsorted || p.out_bam_coord);
```

with:

```rust
    // Decide whether to keep the SAM on disk (`--outSAMtype SAM` or nothing
    // set — default), or to transiently produce SAM and then
    // convert/sort to BAM (`--outSAMtype BAM Unsorted` or
    // `--outSAMtype BAM SortedByCoordinate`). Phase-2a always goes via
    // the SAM intermediate.
    let want_sam = emit_sam && (p.out_sam_bool || (!p.out_bam_unsorted && !p.out_bam_coord));
    let want_bam_unsorted = emit_sam && p.out_bam_unsorted;
    let want_bam_coord = emit_sam && p.out_bam_coord;
    let want_bam = want_bam_unsorted || want_bam_coord;
```

- [ ] **Step 3: Update post-align conversion branch**

Replace the block at lines 394-405 (the `if want_bam { ... }` block):

```rust
    if want_bam_unsorted {
        let compression_level = Some(p.out_bam_compression.clamp(0, 9) as u32);
        star_bam::sam_to_bam::convert_sam_file_to_bam(
            std::path::Path::new(&sam_scratch_path),
            std::path::Path::new(&bam_path),
            compression_level,
        )?;
        let _ = std::fs::remove_file(&sam_scratch_path);
    } else if want_bam_coord {
        let compression_level = Some(p.out_bam_compression.clamp(0, 9) as u32);
        let sorted_bam_path = format!("{out_prefix}Aligned.sortedByCoord.out.bam");
        let tmp_root = if p.out_bam_sort_tmp_dir.is_empty() {
            format!("{out_prefix}_STARtmp_BAMsort")
        } else {
            p.out_bam_sort_tmp_dir.clone()
        };
        std::fs::create_dir_all(&tmp_root)?;

        // Build a per-pass BamOutput, feed the SAM intermediate
        // through it, run the sort, clean up.
        let n_bins = p.out_bam_coord_nbins.max(2);
        let mut bam_out = star_bam::bam_output::BamOutput::new(
            0,
            std::path::Path::new(&tmp_root),
            n_bins,
            p.chunk_out_bam_size_bytes,
        )?;
        star_bam::bam_sort_feed::feed_sam_into_bam_output(
            std::path::Path::new(&sam_scratch_path),
            &mut bam_out,
        )?;
        bam_out.coord_flush()?;

        // Re-read the SAM header bytes to hand to the final writer.
        let sam_header_bytes = std::fs::read(&sam_scratch_path)?;
        // Keep only header lines (prefix lines starting with '@').
        let mut hdr_end = 0usize;
        for line in sam_header_bytes.split_inclusive(|&b| b == b'\n') {
            if !line.starts_with(b"@") { break; }
            hdr_end += line.len();
        }
        let sam_header = &sam_header_bytes[..hdr_end];

        star_bam::bam_sort_coord::sort_bam_by_coordinate(
            std::path::Path::new(&tmp_root),
            /* n_threads */ 1,
            n_bins,
            &[bam_out.bin_total_n().to_vec()],
            &[bam_out.bin_total_bytes().to_vec()],
            sam_header,
            std::path::Path::new(&sorted_bam_path),
            compression_level,
            p.limit_bam_sort_ram,
        )?;
        let _ = std::fs::remove_file(&sam_scratch_path);
        let _ = std::fs::remove_dir_all(&tmp_root);
    }
```

- [ ] **Step 4: Build**

Run: `cargo build -p star-cli 2>&1 | tail -15`
Expected: clean build.

- [ ] **Step 5: Smoke test locally (no reference comparison yet)**

Run from repo root:

```bash
mkdir -p /tmp/star-rs-smoke-sorted
./target/debug/star --runMode alignReads --runThreadN 1 \
  --genomeDir "$(pwd)/tests/fixtures/tiny_genome" \
  --readFilesIn "$(pwd)/tests/fixtures/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /tmp/star-rs-smoke-sorted/
samtools view -c /tmp/star-rs-smoke-sorted/Aligned.sortedByCoord.out.bam
samtools quickcheck /tmp/star-rs-smoke-sorted/Aligned.sortedByCoord.out.bam && echo OK
```

Expected: record count > 0, `quickcheck` prints `OK`. If `tests/fixtures/tiny_genome` is not the right path, look at how `tests/e2e.sh` resolves `$GEN_REF` and reuse that.

- [ ] **Step 6: Commit**

```bash
git add crates/star-cli/src/main.rs
git commit -m "star-cli: wire --outSAMtype BAM SortedByCoordinate via phase-2a pipeline"
```

---

## Task 11: e2e test (single-thread + multi-thread)

**Files:**
- Modify: `tests/e2e.sh` (add section "[9] alignReads --outSAMtype BAM SortedByCoordinate" after line 272).

- [ ] **Step 1: Add mapped/unmapped split comparison helper**

Edit `tests/e2e.sh` — insert after the `cmp_bam_body` function (around line 248):

```bash
cmp_sorted_bam_body() {
  # Phase-2a acceptance: mapped records compared with within-coord qname
  # sort (algorithm parity without iRead coupling); unmapped records
  # compared as a multiset (sort | diff).
  local label="$1" a="$2" b="$3"
  if ! command -v samtools >/dev/null 2>&1; then
    echo "  SKIP $label (samtools not available)"
    return 0
  fi
  local am bm au bu
  am="$(mktemp)"; bm="$(mktemp)"; au="$(mktemp)"; bu="$(mktemp)"
  # Mapped: stable within-coord-bucket order by qname.
  samtools view -F 4 "$a" | awk 'BEGIN{FS=OFS="\t"} {print $3"_"$4, $1, $0}' | sort -k1,1 -k2,2 | cut -f3- > "$am"
  samtools view -F 4 "$b" | awk 'BEGIN{FS=OFS="\t"} {print $3"_"$4, $1, $0}' | sort -k1,1 -k2,2 | cut -f3- > "$bm"
  # Unmapped: multiset.
  samtools view -f 4 "$a" | sort > "$au"
  samtools view -f 4 "$b" | sort > "$bu"
  local ok=1
  if ! cmp -s "$am" "$bm"; then ok=0; echo "  DIFF mapped:"; diff "$am" "$bm" | head -8; fi
  if ! cmp -s "$au" "$bu"; then ok=0; echo "  DIFF unmapped:"; diff "$au" "$bu" | head -8; fi
  if [ "$ok" -eq 1 ]; then pass "$label"; else fail "$label"; fi
  rm -f "$am" "$bm" "$au" "$bu"
}
```

- [ ] **Step 2: Add test block**

Edit `tests/e2e.sh` — insert after line 272 (end of section 8), before the "All regression tests passed." final block:

```bash
echo
echo "[9] alignReads --outSAMtype BAM SortedByCoordinate (phase 2a)"
R9="$SCRATCH/bamsort_rs/"; C9="$SCRATCH/bamsort_ref/"; mkdir -p "$R9" "$C9"
"$STAR_REF" --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$C9" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$R9" >/dev/null 2>&1
cmp_sorted_bam_body "bamsort/Aligned.sortedByCoord.out.bam (semantic)" \
  "$R9/Aligned.sortedByCoord.out.bam" "$C9/Aligned.sortedByCoord.out.bam"

R9m="$SCRATCH/bamsort_rs4/"; C9m="$SCRATCH/bamsort_ref4/"; mkdir -p "$R9m" "$C9m"
"$STAR_REF" --runMode alignReads --runThreadN 4 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$C9m" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 4 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$R9m" >/dev/null 2>&1
cmp_sorted_bam_body "bamsort-mt/Aligned.sortedByCoord.out.bam (semantic)" \
  "$R9m/Aligned.sortedByCoord.out.bam" "$C9m/Aligned.sortedByCoord.out.bam"
```

- [ ] **Step 3: Run e2e tests**

Run (build once, then test):

```bash
cargo build --release -p star-cli 2>&1 | tail -5
bash tests/e2e.sh 2>&1 | tail -40
```

Expected: the two new lines
```
pass bamsort/Aligned.sortedByCoord.out.bam (semantic)
pass bamsort-mt/Aligned.sortedByCoord.out.bam (semantic)
```
are present, and the final line reads `All regression tests passed.`.

If the mapped diff fails, inspect the first few mismatched lines the helper prints — the most likely causes are (a) sort-order within the same `(ref, pos)` bucket producing a CIGAR/QNAME mismatch (fix: refine the sort key in the helper), or (b) missing records (fix: check `BamOutput::coord_flush` is seeing every iRead — add a log line temporarily).

If the unmapped diff fails, the multiset is wrong — usually because the unmapped merge dropped or duplicated records. Add a counter assert in `merge_unmapped_bin` to narrow it down.

- [ ] **Step 4: Commit**

```bash
git add tests/e2e.sh
git commit -m "tests/e2e: add BAM SortedByCoordinate phase-2a semantic comparison"
```

---

## Task 12: Regression-notes + roadmap update

**Files:**
- Modify: `docs/regression-notes.md` (append a note on phase-2a acceptance).
- Modify: `docs/roadmap.md:80` (already updated earlier; verify still accurate).

- [ ] **Step 1: Append phase-2a note**

Append to `docs/regression-notes.md`:

```markdown

## BAM SortedByCoordinate — phase 2a

Implemented as of 2026-04-18. Acceptance is **semantic match via `samtools view`**,
not byte-exact: mapped records are diffed within-`(ref, pos)` sorted by
`qname`; unmapped records are diffed as a multiset. iRead is assigned as
a per-pass monotonic counter (not globally aligned with C++ STAR's
ordinal). The final BAM writer is `noodles-bam`, not the ported
`alignBAM` encoder; byte-exact acceptance waits on T1 / phase 2b.

Verified on `tests/fixtures/tiny/` with `runThreadN=1` and `runThreadN=4`
in `tests/e2e.sh` `[9]`. The sort pipeline itself is single-threaded at
the feeder layer (one `BamOutput`); per-bin sort and unmapped merge are
parallel via rayon.
```

- [ ] **Step 2: Commit**

```bash
git add docs/regression-notes.md
git commit -m "docs: regression-notes for BAM SortedByCoordinate phase 2a"
```

---

## Self-Review

**Spec coverage** (roadmap.md T2 phase-2a bullets):

- [x] `coordOneAlign` — Task 3
- [x] `coordBins` boundary estimation — Task 4
- [x] `coordFlush` — Task 5
- [x] `BAMbinSortByCoordinate` sort key `(coord, iRead, record_offset)` — Task 6 (meta.sort_key)
- [x] `BAMbinSortUnmapped` multi-way merge by iRead — Task 7
- [x] `bamSortByCoordinate` wiring (bin-size aggregation + RAM check + per-bin parallelism + concat) — Task 9
- [x] `noodles-bam` writer (not STAR encoder) — Task 9 `assemble_final_bam`
- [x] Rust-native temp layout — Task 1
- [x] iRead stable-within-thread only — Task 8 (monotonic counter)
- [x] SAM intermediate as input source — Task 8
- [x] `--outBAMsortingThreadN`, `--outBAMsortingBinsN`, `--limitBAMsortRAM` — Task 2
- [x] `ReadAlign_stitchPieces` / `BAMoutput_coordUnmapped` chunk wiring — phase-2a simplification: we feed from the SAM intermediate instead of wiring into the per-chunk alignment path; captured in Task 10 comment "Phase-2a always goes via the SAM intermediate."
- [x] Unit tests `coord_bins_matches_cpp_fixture`, `mapped_bin_sort_key_matches_cpp`, `unmapped_merge_matches_cpp_read_order` — covered by Tasks 4, 6, 7 (tests named slightly differently; same scope).
- [x] Semantic acceptance: mapped bucket-diff + unmapped multiset diff — Task 11

**Type consistency:**
- `BamOutput::new(i_chunk: u32, tmp_root: &Path, n_bins: u32, chunk_bytes: u64)` — used identically in Tasks 3, 5, 9, 10.
- `coord_one_align(&mut self, ref_id: u32, pos: u32, i_read: u64, bam_bytes: Vec<u8>)` — Tasks 3, 8.
- `TmpRecord { ref_id, pos, i_read, bam_bytes }` field order identical across Tasks 1, 6, 7, 8.
- `UNMAPPED_SENTINEL: u32 = u32::MAX` — Tasks 1, 3, 8.
- `sort_bam_by_coordinate` signature — Task 9 + Task 10 call-site match (check when executing Task 10).

**Placeholder scan:** none. Every step has executable commands + complete code.

**Open item flagged in-plan:** Task 8 and Task 9 depend on noodles-bam 0.88 surface details (`bam::Record::try_from(&[u8])`, `write_record`, `try_from_alignment_record`, `Record::as_ref() -> &[u8]`) that may have moved. The tasks explicitly note this and define the contract the agent must preserve if adapting. If the noodles API forces a different serialisation hop, the workaround is to stash the `sam::Record` in `TmpRecord::bam_bytes` (as text bytes), and let the final writer re-parse — slower but always works. Document the choice in the commit message.

---

## Execution Handoff

**Plan complete and saved to `docs/superpowers/plans/2026-04-18-bam-sorted-by-coord-phase-2a.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

**Which approach?**
