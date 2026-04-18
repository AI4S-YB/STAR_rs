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
use std::io::{BufReader, BufWriter, Cursor};
use std::path::Path;

use anyhow::{bail, Context, Result};
use bstr::BString;
use noodles_bam as bam;
use noodles_sam as sam;
use rayon::prelude::*;

use crate::bam_sort_bin::sort_mapped_bin;
use crate::bam_sort_feed::decode_bam_bytes;
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
    //
    // Caveat: phase-2a `sort_mapped_bin` holds the loaded bin bytes AND
    // a parallel `Vec<TmpRecord>` (with cloned bam_bytes) during the
    // sort, so real peak is ~2 * bytes + 24 * n. We under-report here —
    // intentional: mirrors STAR's estimate and keeps `--limitBAMsortRAM`
    // ergonomics consistent with the C++ side. Phase 2b should either
    // tighten this estimate or eliminate the double buffer by sorting
    // indices into the original flat byte buffer.
    let mapped = per_bin_n.len().saturating_sub(1);
    let mut peak = 0u64;
    for ib in 0..mapped {
        let m = per_bin_bytes[ib] + 24 * per_bin_n[ib];
        if m > peak {
            peak = m;
        }
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

    // Parallel fan-out: process the unmapped bin alongside mapped bins in
    // reverse bin order so unmapped (typically largest) starts first.
    let all_bins: Vec<u32> = (0..n_bins_total).rev().collect();
    all_bins.par_iter().try_for_each(|&ibin| -> Result<()> {
        let i = ibin as usize;
        if ibin == (n_bins_total - 1) {
            merge_unmapped_bin(ibin, n_threads, tmp_root)?;
        } else {
            if bin_n[i] == 0 {
                return Ok(());
            }
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
    let mut header = hdr_reader
        .read_header()
        .context("parse SAM header for final BAM")?;

    // Patch SO:coordinate into the header.
    // In noodles-sam 0.84, SortOrder lives in other_fields (not a typed
    // field), accessed via the SORT_ORDER tag constant.
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    if let Some(hd) = header.header_mut().as_mut() {
        hd.other_fields_mut()
            .insert(SORT_ORDER, BString::from("coordinate"));
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
            let bam_rec = decode_bam_bytes(&rec.bam_bytes).with_context(|| {
                format!(
                    "decoding TmpRecord → bam::Record (bin {ibin}, iRead {})",
                    rec.i_read
                )
            })?;
            bam_w.write_record(&header, &bam_rec)?;
        }
        let _ = std::fs::remove_file(&path);
    }
    bam_w.try_finish()?;
    Ok(())
}
