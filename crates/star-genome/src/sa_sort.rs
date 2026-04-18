//! Chunked suffix-array sort driver.
//!
//! Port of `Genome_genomeGenerate.cpp` lines 213-349 (the inner scope after
//! "sort SA chunks"). Algorithm:
//!
//! 1. Flip the genome buffer in place so that for the forward-strand suffix
//!    comparator the "stream" runs backwards from the anchor.
//! 2. Count SA entries bucketed by a 16-bit prefix of 4 consecutive bytes
//!    (`(G[i]<<12) | (G[i-1]<<8) | (G[i-2]<<4) | G[i-3]`).
//! 3. Partition the 65536 buckets into chunks so that each chunk fits in
//!    `sa_chunk_size = (limitGenomeGenerateRAM - nG1alloc) / 8 / runThreadN
//!                     * 6/10`.
//! 4. For each chunk: collect SA indices whose 16-bit prefix falls in the
//!    chunk's bucket range, qsort them using [`compare_suffixes`], and write
//!    to a temp file `SA_<i>`.
//! 5. Read back chunk files in order and pack into the final `SA`
//!    `PackedArray` applying the strand marker:
//!    ```text
//!    packed = (i < nGenome) ? i : ((i - nGenome) | (1 << GstrandBit));
//!    ```
//!
//! This module implements the **single-threaded bit-exact path**. The C++
//! uses OpenMP `#pragma omp parallel for ordered` — threads sort different
//! chunks in parallel but write out in the original chunk order (`ordered`).
//! Because each chunk's sort is deterministic (anti-stable tie-break by
//! index), serial vs. parallel only changes scheduling, not output bytes.
//! Parallelizing with `rayon::par_iter` is the M4 upgrade.
//!
//! # Safety
//! All helpers take the genome buffer as a raw pointer because the
//! comparator reads 7 bytes *before* the anchor (`globalG-7`). Callers must
//! guarantee the 7-byte pre-pad is valid (the C++ relies on positive
//! `G-7` being inside the earlier allocation; we explicitly zero-pad).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Result;
use star_core::packed::PackedArray;

use crate::sa::compare_suffixes;

/// Parameters bag for [`sort_suffix_array`].
pub struct SaSortParams {
    pub n_genome: u64,
    pub sparse_d: u64,
    pub suffix_length_max: u64,
    pub limit_genome_generate_ram: u64,
    pub n_g1_alloc: u64,
    pub run_thread_n: u32,
}

/// Stage 1: Count 16-bit-prefix buckets. Operates on the already-flipped
/// genome (caller must do `G[2N-1-i] <-> G[i]` swap per line 216).
///
/// `g_ptr` must point at the logical `G[0]` position and be backed by ≥3 bytes
/// of pre-padding (so `G[-3..0]` reads are valid, matching the C++
/// allocation).
///
/// # Safety
/// Caller guarantees `g_ptr` is valid for reads of `[-3, 2*n_genome)`.
pub unsafe fn count_prefix_buckets(
    g_ptr: *const u8,
    n_genome: u64,
    sparse_d: u64,
) -> (Vec<u64>, u64) {
    const IND_PREF_N: usize = 1 << 16;
    let mut ind_pref_count = vec![0u64; IND_PREF_N];
    let mut n_sa: u64 = 0;
    let mut ii: u64 = 0;
    while ii < 2 * n_genome {
        unsafe {
            let cur = *g_ptr.add(ii as usize);
            if cur < 4 {
                let p1 = ((cur as usize) << 12)
                    + ((*g_ptr.offset(ii as isize - 1) as usize) << 8)
                    + ((*g_ptr.offset(ii as isize - 2) as usize) << 4)
                    + (*g_ptr.offset(ii as isize - 3) as usize);
                ind_pref_count[p1] += 1;
                n_sa += 1;
            }
        }
        ii += sparse_d;
    }
    (ind_pref_count, n_sa)
}

/// Stage 2: Partition buckets into chunks of at most `sa_chunk_size` entries.
/// Returns `(ind_pref_start, ind_pref_chunk_count)` matching the C++ layout.
pub fn partition_chunks(
    ind_pref_count: &[u64],
    n_sa: u64,
    sa_chunk_size: u64,
) -> (Vec<u64>, Vec<u64>) {
    let cap = n_sa / sa_chunk_size + 1;
    let mut ind_pref_start = Vec::with_capacity((cap * 2) as usize);
    let mut ind_pref_chunk_count = Vec::with_capacity((cap * 2) as usize);
    ind_pref_start.push(0u64);
    let mut chunk_size1: u64 = ind_pref_count[0];
    for ii in 1..ind_pref_count.len() {
        chunk_size1 += ind_pref_count[ii];
        if chunk_size1 > sa_chunk_size {
            ind_pref_start.push(ii as u64);
            ind_pref_chunk_count.push(chunk_size1 - ind_pref_count[ii]);
            chunk_size1 = ind_pref_count[ii];
        }
    }
    ind_pref_chunk_count.push(chunk_size1);
    ind_pref_start.push(ind_pref_count.len() as u64 + 1);
    (ind_pref_start, ind_pref_chunk_count)
}

/// Stage 3: Sort a single chunk (no I/O). Collects SA indices whose prefix
/// falls in `[bucket_lo, bucket_hi)`, sorts via [`compare_suffixes`], and
/// returns the sorted index list.
///
/// # Safety
/// See [`count_prefix_buckets`].
pub unsafe fn sort_chunk(
    g_ptr: *const u8,
    n_genome: u64,
    sparse_d: u64,
    suffix_length_max: u64,
    bucket_lo: u64,
    bucket_hi: u64,
    chunk_count: u64,
) -> Vec<u64> {
    let mut out = Vec::with_capacity(chunk_count as usize);
    let mut ii: u64 = 0;
    while ii < 2 * n_genome {
        unsafe {
            let cur = *g_ptr.add(ii as usize);
            if cur < 4 {
                let p1 = ((cur as u64) << 12)
                    + ((*g_ptr.offset(ii as isize - 1) as u64) << 8)
                    + ((*g_ptr.offset(ii as isize - 2) as u64) << 4)
                    + (*g_ptr.offset(ii as isize - 3) as u64);
                if p1 >= bucket_lo && p1 < bucket_hi {
                    out.push(ii);
                }
            }
        }
        ii += sparse_d;
    }
    debug_assert_eq!(out.len() as u64, chunk_count);

    let global_l = suffix_length_max / 8;
    // Sort via unsafe comparator. We intentionally use `sort_by` (stable) —
    // the comparator itself is anti-stable by SA index, so stable vs unstable
    // wrapping does not affect final order.
    out.sort_by(|a, b| unsafe { compare_suffixes(g_ptr, global_l, *a, *b) });
    // Apply `saChunk[ii] = 2*nGenome - 1 - saChunk[ii]`.
    for v in &mut out {
        *v = 2 * n_genome - 1 - *v;
    }
    out
}

/// Stage 4: Write a sorted chunk to disk in raw `u64` little-endian form.
pub fn write_chunk(dir: &Path, idx: u64, chunk: &[u64]) -> Result<PathBuf> {
    let path = dir.join(format!("SA_{idx}"));
    let mut f = BufWriter::new(File::create(&path)?);
    // 1:1 with C++ fstreamWriteBig(saChunk, sizeof(uint)*count) where
    // uint=u64 on x86-64.
    for &v in chunk {
        f.write_all(&v.to_le_bytes())?;
    }
    f.flush()?;
    Ok(path)
}

/// Stage 5: Pack chunks into final SA PackedArray.
pub fn pack_chunks(
    chunk_paths: &[PathBuf],
    sa: &mut PackedArray,
    n_genome: u64,
    g_strand_bit: u8,
) -> Result<u64> {
    let n2_bit: u64 = 1 << g_strand_bit;
    let mut packed_ind: u64 = 0;
    let mut buf = vec![0u8; 10_000_000 * 8];
    for path in chunk_paths {
        let mut f = File::open(path)?;
        loop {
            use std::io::Read;
            let n = f.read(&mut buf)?;
            if n == 0 {
                break;
            }
            let n_elems = n / 8;
            for i in 0..n_elems {
                let bytes = &buf[i * 8..i * 8 + 8];
                let v = u64::from_le_bytes(bytes.try_into().unwrap());
                let packed = if v < n_genome {
                    v
                } else {
                    (v - n_genome) | n2_bit
                };
                sa.write_packed(packed_ind + i as u64, packed);
            }
            packed_ind += n_elems as u64;
        }
        let _ = std::fs::remove_file(path);
    }
    Ok(packed_ind)
}

/// Top-level: orchestrate all stages. Deletes `SA_<i>` chunk files on the
/// way out. Returns `nSA` for caller sanity-checking.
///
/// # Safety
/// `g_ptr` must point at the logical `G[0]` of the flipped genome with
/// ≥3 bytes of pre-padding and ≥8 bytes of post-padding. See
/// [`count_prefix_buckets`].
pub unsafe fn sort_suffix_array(
    g_ptr: *const u8,
    sa: &mut PackedArray,
    params: &SaSortParams,
    g_strand_bit: u8,
    dir: &Path,
) -> Result<u64> {
    let (counts, n_sa) = unsafe { count_prefix_buckets(g_ptr, params.n_genome, params.sparse_d) };

    let sa_chunk_size = {
        let raw =
            (params.limit_genome_generate_ram - params.n_g1_alloc) / 8 / params.run_thread_n as u64;
        let raw = raw * 6 / 10;
        if params.run_thread_n > 1 {
            raw.min(n_sa / (params.run_thread_n as u64 - 1))
        } else {
            raw
        }
    };

    let (ind_pref_start, chunk_counts) = partition_chunks(&counts, n_sa, sa_chunk_size);

    let mut paths = Vec::with_capacity(chunk_counts.len());
    for (i, &cnt) in chunk_counts.iter().enumerate() {
        let sorted = unsafe {
            sort_chunk(
                g_ptr,
                params.n_genome,
                params.sparse_d,
                params.suffix_length_max,
                ind_pref_start[i],
                ind_pref_start[i + 1],
                cnt,
            )
        };
        paths.push(write_chunk(dir, i as u64, &sorted)?);
    }

    let packed = pack_chunks(&paths, sa, params.n_genome, g_strand_bit)?;
    if packed != n_sa {
        anyhow::bail!("SA pack mismatch: packed={packed} expected={n_sa}");
    }
    Ok(n_sa)
}
