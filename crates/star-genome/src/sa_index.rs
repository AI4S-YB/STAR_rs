//! 1:1 port of `genomeSAindex.cpp`. Builds the prefix-truncated SAindex
//! given the already-sorted suffix array.

use anyhow::Result;
use star_core::packed::PackedArray;

use crate::genome::Genome;
use crate::sa_funs::fun_calc_sai_from_sa;

/// Port of `genomeSAindex()` (genomeSAindex.cpp:6).
///
/// Initializes `Genome::{genomeSAindexStart,nSAi,SAiMark*,GstrandMask}` then
/// builds `SAi` PackedArray. Caller owns `SA` (already built & packed by
/// `sa_sort::sort_suffix_array`).
///
/// # Safety
/// `g_ptr` must point at the logical `G[0]` with ≥16 bytes of pre-padding.
pub unsafe fn genome_sa_index(
    g_ptr: *const u8,
    sa: &PackedArray,
    sa_index_nbases: u32,
    map_gen: &mut Genome,
) -> Result<PackedArray> {
    map_gen.genome_sa_index_start = Vec::with_capacity(sa_index_nbases as usize + 1);
    map_gen.genome_sa_index_start.push(0);
    for ii in 1..=sa_index_nbases as u64 {
        let prev = *map_gen.genome_sa_index_start.last().unwrap();
        map_gen.genome_sa_index_start.push(prev + (1u64 << (2 * ii)));
    }
    map_gen.n_sai = map_gen.genome_sa_index_start[sa_index_nbases as usize];

    map_gen.sai_mark_n_bit = map_gen.g_strand_bit + 1;
    map_gen.sai_mark_absent_bit = map_gen.g_strand_bit + 2;
    map_gen.sai_mark_n_mask_c = 1u64 << map_gen.sai_mark_n_bit;
    map_gen.sai_mark_n_mask = !map_gen.sai_mark_n_mask_c;
    map_gen.sai_mark_absent_mask_c = 1u64 << map_gen.sai_mark_absent_bit;
    map_gen.sai_mark_absent_mask = !map_gen.sai_mark_absent_mask_c;

    let mut sai = PackedArray::default();
    sai.define_bits((map_gen.g_strand_bit + 3) as u32, map_gen.n_sai);
    sai.allocate_array();

    unsafe {
        genome_sa_index_chunk(g_ptr, sa, &mut sai, 0, sa.length.saturating_sub(1),
            sa_index_nbases, map_gen)?;
    }
    Ok(sai)
}

/// Port of `genomeSAindexChunk()` (genomeSAindex.cpp:117).
unsafe fn genome_sa_index_chunk(
    g_ptr: *const u8,
    sa: &PackedArray,
    sai: &mut PackedArray,
    i_sa1: u64,
    i_sa2: u64,
    sa_index_nbases: u32,
    map_gen: &Genome,
) -> Result<()> {
    // `ind0[iL] = -1` (unsigned overflow) per C++.
    let mut ind0 = vec![u64::MAX; sa_index_nbases as usize];

    let isa_step = map_gen.n_sa / (1u64 << (2 * sa_index_nbases)) + 1;

    let mut isa = i_sa1;
    let (mut ind_full, mut il4) = unsafe {
        fun_calc_sai_from_sa(g_ptr, sa, map_gen, isa, sa_index_nbases as i32)
    };

    while isa <= i_sa2 {
        for i_l in 0..sa_index_nbases as i32 {
            let ind_pref = ind_full >> (2 * (sa_index_nbases as i32 - 1 - i_l));

            if i_l == il4 {
                // Suffix contains N: mark the last seen prefix with N-mask.
                for i_l1 in i_l..sa_index_nbases as i32 {
                    let start = map_gen.genome_sa_index_start[i_l1 as usize];
                    let idx0 = ind0[i_l1 as usize];
                    let v = sai.get(start + idx0) | map_gen.sai_mark_n_mask_c;
                    sai.write_packed(start + idx0, v);
                }
                break;
            }

            let start = map_gen.genome_sa_index_start[i_l as usize];
            let prev = ind0[i_l as usize];

            if ind_pref > prev || isa == 0 {
                sai.write_packed(start + ind_pref, isa);
                let lo = if prev == u64::MAX { 0 } else { prev + 1 };
                for ii in lo..ind_pref {
                    sai.write_packed(start + ii, isa | map_gen.sai_mark_absent_mask_c);
                }
                ind0[i_l as usize] = ind_pref;
            } else if ind_pref < prev {
                anyhow::bail!(
                    "BUG: next index is smaller than previous at iL={i_l}, isa={isa}"
                );
            }
        }

        // Advance isa with large-step + binary search.
        unsafe {
            fun_sai_find_next_index(g_ptr, sa, isa_step, &mut isa, &mut ind_full,
                &mut il4, sa_index_nbases as i32, map_gen);
        }
    }

    // Fill trailing empty slots.
    for i_l in 0..sa_index_nbases as usize {
        let start = map_gen.genome_sa_index_start[i_l];
        let end = map_gen.genome_sa_index_start[i_l + 1];
        let prev = ind0[i_l];
        let lo = if prev == u64::MAX { start } else { start + prev + 1 };
        for ii in lo..end {
            sai.write_packed(ii, map_gen.n_sa | map_gen.sai_mark_absent_mask_c);
        }
    }
    Ok(())
}

/// Port of `funSAiFindNextIndex()` (genomeSAindex.cpp:178).
unsafe fn fun_sai_find_next_index(
    g_ptr: *const u8,
    sa: &PackedArray,
    isa_step: u64,
    isa: &mut u64,
    ind_full: &mut u64,
    il4: &mut i32,
    l: i32,
    map_gen: &Genome,
) {
    let ind_full_prev = *ind_full;
    let il4_prev = *il4;
    *isa += isa_step;
    while *isa < map_gen.n_sa {
        let (nf, nil) = unsafe { fun_calc_sai_from_sa(g_ptr, sa, map_gen, *isa, l) };
        *ind_full = nf;
        *il4 = nil;
        if !(nf == ind_full_prev && nil == il4_prev) {
            break;
        }
        *isa += isa_step;
    }

    if *isa >= map_gen.n_sa {
        let (nf, nil) = unsafe {
            fun_calc_sai_from_sa(g_ptr, sa, map_gen, map_gen.n_sa - 1, l)
        };
        *ind_full = nf;
        *il4 = nil;
        if nf == ind_full_prev && nil == il4_prev {
            *isa = map_gen.n_sa;
            return;
        }
    }

    let mut i1 = isa.saturating_sub(isa_step);
    let mut i2 = (*isa).min(map_gen.n_sa - 1);
    while i1 + 1 < i2 {
        *isa = i1 / 2 + i2 / 2 + (i1 % 2 + i2 % 2) / 2;
        let (nf, nil) = unsafe { fun_calc_sai_from_sa(g_ptr, sa, map_gen, *isa, l) };
        *ind_full = nf;
        *il4 = nil;
        if nf == ind_full_prev && nil == il4_prev {
            i1 = *isa;
        } else {
            i2 = *isa;
        }
    }
    if *isa == i1 {
        *isa = i2;
        let (nf, nil) = unsafe { fun_calc_sai_from_sa(g_ptr, sa, map_gen, *isa, l) };
        *ind_full = nf;
        *il4 = nil;
    }
}
