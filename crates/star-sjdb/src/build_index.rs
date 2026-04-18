//! 1:1 port of `sjdbBuildIndex.cpp` + `Genome::insertSequences`.
//!
//! Inserts the junction sequences (`Gsj`, produced by `sjdb_prepare`) into
//! the running genome's G / SA / SAi. Output: an updated `Genome` whose
//! `sjdb_n`, `n_genome`, `n_sa`, `n_sa_byte`, `sj_g_start` have been
//! advanced, and SA/SAi reflect the inserted junctions.
//!
//! The incremental build keeps the tie-break ordering used by C++ so that
//! pass-2 alignments match bit-exact.

use star_core::packed::PackedArray;
use star_core::seq::complement_seq_numbers;
use star_core::types::GENOME_SPACING_CHAR;
use star_genome::genome::Genome;
use star_genome::load::LOAD_L;
use star_genome::sa_funs::{fun_calc_sai, fun_compare_uint_and_suffixes, suffix_array_search1};
use star_params::parameters::Parameters;

use crate::binary_search2::binary_search2_pair;

/// Port of `sjdbBuildIndex` (sjdbBuildIndex.cpp:16-333).
///
/// - `map_gen` is the genome being extended; on entry, its `sjdb_*` fields
///   already reflect the combined sjdb from `sjdb_prepare`.
/// - `map_gen1` is the pre-insertion snapshot used to locate "old" junctions.
/// - `gsj` is the Gsj buffer produced by `sjdb_prepare` (length
///   `2 * sjdb_length * sjdb_n + 1`).
/// - Caller must have resized `map_gen.g` to accommodate the new
///   `chr_start[n_chr_real] + n_gsj` bytes; this function fills the tail.
pub fn sjdb_build_index(
    p: &Parameters,
    gsj: &mut [u8],
    map_gen: &mut Genome,
    map_gen1: &Genome,
) -> anyhow::Result<()> {
    if map_gen.sjdb_n == 0 {
        return Ok(());
    }

    let sjdb_length = map_gen.sjdb_length;
    let sjdb_n = map_gen.sjdb_n;
    let n_gsj = (sjdb_length * sjdb_n) as usize;

    // Close each junction with a spacer char, then RC-mirror into the 2nd half.
    for ii in 1..=sjdb_n {
        gsj[(ii * sjdb_length) as usize - 1] = GENOME_SPACING_CHAR;
    }
    gsj[n_gsj * 2] = GENOME_SPACING_CHAR;
    for ii in 0..n_gsj {
        let v = gsj[ii];
        gsj[n_gsj * 2 - 1 - ii] = if v < 4 { 3 - v } else { v };
    }

    // G1c: complement (not reverse) of Gsj, used as `s2[1]` in
    // suffix_array_search1. +1 for the end spacer.
    let mut g1c = vec![0u8; n_gsj * 2 + 1];
    complement_seq_numbers(&gsj[..n_gsj * 2 + 1], &mut g1c);

    // oldSJind[old_sjdb_index] = merged position in new sjdb_*.
    let mut old_sj_ind = vec![0u32; map_gen1.sjdb_n as usize];

    let n_indices_sj1 = sjdb_length;
    let total_pairs = (2 * sjdb_n * n_indices_sj1) as usize;
    let mut ind_array: Vec<[u64; 2]> = vec![[0u64; 2]; total_pairs + 1];

    let mut sj_new_count: u64 = 0;
    let chr_start_end = map_gen.chr_start[map_gen.n_chr_real as usize];

    for isj in 0..2 * sjdb_n {
        let isj_us = isj as usize;
        let seq0 = &gsj[isj_us * sjdb_length as usize..];
        let seq1 = &g1c[isj_us * sjdb_length as usize..];

        let isj1 = if isj < sjdb_n {
            isj
        } else {
            2 * sjdb_n - 1 - isj
        };
        let sjdb_ind = if map_gen1.sjdb_n == 0 {
            -1i64
        } else {
            binary_search2_pair(
                map_gen.sjdb_start[isj1 as usize],
                map_gen.sjdb_end[isj1 as usize],
                &map_gen1.sjdb_start,
                &map_gen1.sjdb_end,
            )
        };
        if sjdb_ind < 0 {
            sj_new_count += 1;
        } else {
            old_sj_ind[sjdb_ind as usize] = isj1 as u32;
        }

        for istart1 in 0..n_indices_sj1 {
            let istart = istart1;
            let ind1 = (isj * n_indices_sj1 + istart1) as usize;
            if sjdb_ind >= 0 || seq0[istart as usize] > 3 {
                ind_array[ind1][0] = u64::MAX;
            } else {
                unsafe {
                    ind_array[ind1][0] = suffix_array_search1(
                        map_gen,
                        &map_gen.sa,
                        [seq0, seq1],
                        istart,
                        10000,
                        u64::MAX,
                        true,
                        0,
                        map_gen.n_sa - 1,
                        0,
                    );
                }
                ind_array[ind1][1] = isj * sjdb_length + istart;
            }
        }
    }
    let sj_new = sj_new_count / 2;

    // Remove invalid entries.
    let mut n_ind = 0usize;
    for ii in 0..total_pairs {
        if ind_array[ii][0] != u64::MAX {
            ind_array[n_ind] = ind_array[ii];
            n_ind += 1;
        }
    }

    // Stable sort with tie-break on suffix ordering.
    let gsj_for_cmp = gsj.to_vec();
    ind_array[..n_ind].sort_by(|a, b| fun_compare_uint_and_suffixes(&gsj_for_cmp, a, b));

    // Sentinel.
    ind_array[n_ind][0] = u64::MAX - 999;
    ind_array[n_ind][1] = u64::MAX - 999;

    // Resize G, SA; track new counts.
    let new_n_genome = chr_start_end + n_gsj as u64;
    map_gen.n_genome = new_n_genome;
    map_gen.n_sa += n_ind as u64;

    let g_strand_bit1 = ((map_gen.n_genome as f64).log2().floor() as u64 + 1).max(32) as u8;
    if g_strand_bit1 > map_gen.g_strand_bit {
        anyhow::bail!(
            "EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem"
        );
    }

    // Allocate SA2 (new suffix array) and fill it by merging old SA with
    // new entries in ind_array.
    let mut sa2 = PackedArray::new();
    sa2.define_bits((map_gen.g_strand_bit + 1) as u32, map_gen.n_sa);
    sa2.allocate_array();

    let n_gsj_new = sj_new * sjdb_length;
    let n2_bit: u64 = 1u64 << map_gen.g_strand_bit;
    let strand_mask: u64 = !n2_bit;

    let mut isj = 0usize;
    let mut isa2: u64 = 0;
    for isa in 0..map_gen1.n_sa {
        while ind_array[isj][0] == isa {
            let mut ind1 = ind_array[isj][1];
            if ind1 < n_gsj as u64 {
                ind1 += chr_start_end;
            } else {
                ind1 = (ind1 - n_gsj as u64) | n2_bit;
            }
            sa2.write_packed(isa2, ind1);
            isa2 += 1;
            isj += 1;
        }

        let mut ind1 = map_gen.sa.get(isa);
        if (ind1 & n2_bit) > 0 {
            let mut ind1s = map_gen1.n_genome - (ind1 & strand_mask);
            if ind1s >= chr_start_end {
                let sj1 = (ind1s - chr_start_end) / sjdb_length;
                let delta = old_sj_ind[sj1 as usize] as i64 - sj1 as i64;
                ind1s = (ind1s as i64 + delta * sjdb_length as i64) as u64;
                ind1 = (map_gen.n_genome - ind1s) | n2_bit;
            } else {
                ind1 += n_gsj_new;
            }
        } else if ind1 >= chr_start_end {
            let sj1 = (ind1 - chr_start_end) / sjdb_length;
            let delta = old_sj_ind[sj1 as usize] as i64 - sj1 as i64;
            ind1 = (ind1 as i64 + delta * sjdb_length as i64) as u64;
        }
        sa2.write_packed(isa2, ind1);
        isa2 += 1;
    }
    while isj < n_ind {
        let mut ind1 = ind_array[isj][1];
        if ind1 < n_gsj as u64 {
            ind1 += chr_start_end;
        } else {
            ind1 = (ind1 - n_gsj as u64) | n2_bit;
        }
        sa2.write_packed(isa2, ind1);
        isa2 += 1;
        isj += 1;
    }

    // ---------- SAi rewrite ----------
    let sa_index_nbases = p.p_ge.g_sa_index_nbases as usize;
    for i_l in 0..sa_index_nbases {
        let mut i_sj = 0usize;
        let start = map_gen.genome_sa_index_start[i_l];
        let end = map_gen.genome_sa_index_start[i_l + 1];
        let mut ind0 = start.wrapping_sub(1);
        for ii in start..end {
            let i_sa1 = map_gen.sai.get(ii);
            let i_sa2 = i_sa1 & map_gen.sai_mark_n_mask & map_gen.sai_mark_absent_mask;

            if i_sj < n_ind && (i_sa1 & map_gen.sai_mark_absent_mask_c) > 0 {
                // Absent index in the old SA — scan for a new junction
                // whose SAi prefix equals `(ii - start)`.
                let i_sj1 = i_sj;
                let gsj_start = ind_array[i_sj][1] as usize;
                let mut ind1 = fun_calc_sai(&gsj[gsj_start..], (i_l + 1) as u32);
                let target = (ii - start) as i64;
                while ind1 < target && (ind_array[i_sj][0] as i64 - 1) < i_sa2 as i64 {
                    i_sj += 1;
                    let gsj_start = ind_array[i_sj][1] as usize;
                    ind1 = fun_calc_sai(&gsj[gsj_start..], (i_l + 1) as u32);
                }
                if ind1 == target {
                    map_gen
                        .sai
                        .write_packed(ii, ind_array[i_sj][0] - 1 + i_sj as u64 + 1);
                    for ii0 in ind0.wrapping_add(1)..ii {
                        map_gen.sai.write_packed(
                            ii0,
                            (ind_array[i_sj][0] - 1 + i_sj as u64 + 1)
                                | map_gen.sai_mark_absent_mask_c,
                        );
                    }
                    i_sj += 1;
                    ind0 = ii;
                } else {
                    i_sj = i_sj1;
                }
            } else {
                while i_sj < n_ind && (ind_array[i_sj][0] as i64 - 1 + 1) < i_sa2 as i64 {
                    i_sj += 1;
                }
                while i_sj < n_ind && (ind_array[i_sj][0] as i64 - 1 + 1) == i_sa2 as i64 {
                    let gsj_start = ind_array[i_sj][1] as usize;
                    let calc = fun_calc_sai(&gsj[gsj_start..], (i_l + 1) as u32);
                    if calc >= (ii - start) as i64 {
                        break;
                    }
                    i_sj += 1;
                }
                map_gen.sai.write_packed(ii, i_sa1 + i_sj as u64);
                for ii0 in ind0.wrapping_add(1)..ii {
                    map_gen
                        .sai
                        .write_packed(ii0, (i_sa2 + i_sj as u64) | map_gen.sai_mark_absent_mask_c);
                }
                ind0 = ii;
            }
        }
    }

    // Mark N-containing prefixes in SAi.
    let sai_mark_n_mask_c = !map_gen.sai_mark_n_mask;
    for isj in 0..n_ind {
        let mut ind1: i64 = 0;
        let gsj_start = ind_array[isj][1] as usize;
        for i_l in 0..sa_index_nbases {
            let g = gsj[gsj_start + i_l] as u64;
            ind1 <<= 2;
            if g > 3 {
                for i_l1 in i_l..sa_index_nbases {
                    ind1 += 3;
                    let mut ind2 = map_gen.genome_sa_index_start[i_l1] as i64 + ind1;
                    while ind2 >= 0
                        && (map_gen.sai.get(ind2 as u64) & map_gen.sai_mark_absent_mask_c) != 0
                    {
                        ind2 -= 1;
                    }
                    if ind2 >= 0 {
                        let cur = map_gen.sai.get(ind2 as u64);
                        map_gen
                            .sai
                            .write_packed(ind2 as u64, cur | sai_mark_n_mask_c);
                    }
                    ind1 <<= 2;
                }
                break;
            } else {
                ind1 += g as i64;
            }
        }
    }

    // Commit SA2 as the new SA.
    map_gen.sa = sa2;
    map_gen.n_sa_byte = map_gen.sa.length_byte;
    map_gen.sj_g_start = chr_start_end;

    // Insert Gsj into G at chr_start[n_chr_real]; leave G1 padding intact.
    let dst = LOAD_L + chr_start_end as usize;
    if dst + n_gsj > map_gen.g.len() {
        // Resize g to accommodate new n_genome + padding.
        let required = dst + n_gsj + LOAD_L;
        if required > map_gen.g.len() {
            map_gen.g.resize(required, GENOME_SPACING_CHAR);
        }
    }
    map_gen.g[dst..dst + n_gsj].copy_from_slice(&gsj[..n_gsj]);

    Ok(())
}
