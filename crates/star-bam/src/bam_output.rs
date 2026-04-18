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
    /// picks the largest index `ib` such that `bin_starts[ib] <= align_g`.
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
