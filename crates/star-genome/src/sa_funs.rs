//! 1:1 port of `SuffixArrayFuns.cpp`.

use crate::genome::{Genome, G_OFFSET};
use star_core::packed::PackedArray;
use star_core::types::GENOME_SPACING_CHAR;

/// Port of `funCalcSAiFromSA` (SuffixArrayFuns.cpp:353).
///
/// Computes the `L`-nucleotide prefix of the suffix at `SA[isa]`. Returns
/// the 2-bit-per-nt packed prefix value and sets `il4` to the first index
/// containing an `N` (or `-1` if none).
///
/// # Safety
/// `g_ptr` must point at the logical `G[0]` with ≥16 bytes of pre-padding
/// (matching `G1 + 100` from C++ `genomeSequenceAllocate`).
pub unsafe fn fun_calc_sai_from_sa(
    g_ptr: *const u8,
    g_sa: &PackedArray,
    map_gen: &Genome,
    isa: u64,
    l: i32,
) -> (u64, i32) {
    let sa_str_raw = g_sa.get(isa);
    let dir_g = (sa_str_raw >> map_gen.g_strand_bit) == 0;
    let sa_str = sa_str_raw & map_gen.g_strand_mask;

    let mut il4: i32 = -1;
    let mut saind: u64 = 0;

    if dir_g {
        for ii in 0..l {
            let off = sa_str as isize + ii as isize;
            let g2 = unsafe { *g_ptr.offset(off) };
            if g2 > 3 {
                il4 = ii;
                saind <<= 2 * (l - ii) as u32;
                return (saind, il4);
            }
            saind = (saind << 2) + g2 as u64;
        }
    } else {
        let base = map_gen.n_genome as isize - sa_str as isize - 1;
        for ii in 0..l {
            let off = base - ii as isize;
            let g2 = unsafe { *g_ptr.offset(off) };
            if g2 > 3 {
                il4 = ii;
                saind <<= 2 * (l - ii) as u32;
                return (saind, il4);
            }
            saind = (saind << 2) + (3 - g2 as u64);
        }
    }
    (saind, il4)
}

/// Median of two unsigned ints (`(a+b)/2` without overflow).
/// Port of `medianUint2` (SuffixArrayFuns.cpp:4).
#[inline]
pub fn median_uint2(a: u64, b: u64) -> u64 {
    a / 2 + b / 2 + (a % 2 + b % 2) / 2
}

/// Port of `compareSeqToGenome` (SuffixArrayFuns.cpp:10).
///
/// Given suffix `SA[isa]`, compare the `N`-long read piece starting at `S`
/// (forward or reverse) against the genome and return the number of matching
/// bases. `comp_res` is set `true` iff the read > genome at the mismatch.
///
/// # Safety
/// `g_ptr` must point at `G[0]` with sufficient padding (matches
/// `LOAD_L = 200`). `read` must provide `read[0]` (forward nt[0..4]) and
/// `read[1]` (rev-complement nt[0..4]) of length ≥ `S + N` (if `dir_r`) or
/// ≥ `S` (if not).
pub unsafe fn compare_seq_to_genome(
    map_gen: &Genome,
    read: &[&[u8]; 2],
    s: u64,
    n: u64,
    l: u64,
    i_sa: u64,
    dir_r: bool,
) -> (u64, bool) {
    let g_ptr = map_gen.g_ptr();
    let sa_str_raw = map_gen.sa.get(i_sa);
    let dir_g = (sa_str_raw >> map_gen.g_strand_bit) == 0;
    let sa_str = sa_str_raw & map_gen.g_strand_mask;

    let limit = (n - l) as isize;

    if dir_r && dir_g {
        let s_base = read[0].as_ptr().add((s + l) as usize);
        let g_base = g_ptr.offset(sa_str as isize + l as isize);
        for ii in 0..limit {
            let sv = *s_base.offset(ii);
            let gv = *g_base.offset(ii);
            if sv != gv {
                return (ii as u64 + l, sv > gv);
            }
        }
        (n, false)
    } else if dir_r && !dir_g {
        let s_base = read[1].as_ptr().add((s + l) as usize);
        let g_base = g_ptr.offset(map_gen.n_genome as isize - 1 - sa_str as isize - l as isize);
        for ii in 0..limit {
            let sv = *s_base.offset(ii);
            let gv = *g_base.offset(-ii);
            if sv != gv {
                let comp = !(sv > gv || gv > 3);
                return (ii as u64 + l, comp);
            }
        }
        (n, false)
    } else if !dir_r && dir_g {
        let s_base = read[1].as_ptr().add((s as isize - l as isize) as usize);
        let g_base = g_ptr.offset(sa_str as isize + l as isize);
        for ii in 0..limit {
            let sv = *s_base.offset(-ii);
            let gv = *g_base.offset(ii);
            if sv != gv {
                return (ii as u64 + l, sv > gv);
            }
        }
        (n, false)
    } else {
        let s_base = read[0].as_ptr().add((s as isize - l as isize) as usize);
        let g_base = g_ptr.offset(map_gen.n_genome as isize - 1 - sa_str as isize - l as isize);
        for ii in 0..limit {
            let sv = *s_base.offset(-ii);
            let gv = *g_base.offset(-ii);
            if sv != gv {
                let comp = !(sv > gv || gv > 3);
                return (ii as u64 + l, comp);
            }
        }
        (n, false)
    }
}

/// Port of `findMultRange` (SuffixArrayFuns.cpp:106).
///
/// Walks the SA range to find the farthest index with identity length `L3`.
#[allow(clippy::too_many_arguments)]
pub unsafe fn find_mult_range(
    map_gen: &Genome,
    i3: u64,
    l3: u64,
    mut i1: u64,
    l1: u64,
    mut i1a: u64,
    l1a: u64,
    mut i1b: u64,
    mut l1b: u64,
    read: &[&[u8]; 2],
    dir_r: bool,
    s: u64,
) -> u64 {
    if l1 < l3 {
        l1b = l1;
        i1b = i1;
        i1a = i3;
    } else if l1a < l1 {
        l1b = l1a;
        i1b = i1a;
        i1a = i1;
    }
    while (i1b + 1 < i1a) || (i1b > i1a + 1) {
        let i1c = median_uint2(i1a, i1b);
        let (l1c, _) = compare_seq_to_genome(map_gen, read, s, l3, l1b, i1c, dir_r);
        if l1c == l3 {
            i1a = i1c;
        } else {
            i1b = i1c;
            l1b = l1c;
        }
        let _ = i1;
        i1 = i1b; // keep i1 in-scope (unused after loop)
    }
    i1a
}

/// Port of `maxMappableLength` (SuffixArrayFuns.cpp:133).
///
/// Binary search in SA space; returns (nrep, updated_L, [i1, i2]).
pub unsafe fn max_mappable_length(
    map_gen: &Genome,
    read: &[&[u8]; 2],
    s: u64,
    n: u64,
    mut i1: u64,
    mut i2: u64,
    dir_r: bool,
    l_in: u64,
) -> (u64, u64, [u64; 2]) {
    let (mut l1, _) = compare_seq_to_genome(map_gen, read, s, n, l_in, i1, dir_r);
    let (mut l2, _) = compare_seq_to_genome(map_gen, read, s, n, l_in, i2, dir_r);
    let mut l = std::cmp::min(l1, l2);

    let (mut l1a, mut i1a) = (l1, i1);
    let (mut l1b, mut i1b) = (l1, i1);
    let (mut l2a, mut i2a) = (l2, i2);
    let (mut l2b, mut i2b) = (l2, i2);

    let mut i3 = i1;
    let mut l3 = l1;
    while i1 + 1 < i2 {
        i3 = median_uint2(i1, i2);
        let (l3c, comp_res) = compare_seq_to_genome(map_gen, read, s, n, l, i3, dir_r);
        l3 = l3c;
        if l3 == n {
            break;
        }
        if comp_res {
            if l3 > l1 {
                l1b = l1a;
                l1a = l1;
                i1b = i1a;
                i1a = i1;
            }
            i1 = i3;
            l1 = l3;
        } else {
            if l3 > l2 {
                l2b = l2a;
                l2a = l2;
                i2b = i2a;
                i2a = i2;
            }
            i2 = i3;
            l2 = l3;
        }
        l = std::cmp::min(l1, l2);
    }

    if l3 < n {
        if l1 > l2 {
            i3 = i1;
            l3 = l1;
        } else {
            i3 = i2;
            l3 = l2;
        }
    }

    let i1_out = find_mult_range(map_gen, i3, l3, i1, l1, i1a, l1a, i1b, l1b, read, dir_r, s);
    let i2_out = find_mult_range(map_gen, i3, l3, i2, l2, i2a, l2a, i2b, l2b, read, dir_r, s);

    (i2_out - i1_out + 1, l3, [i1_out, i2_out])
}

/// Port of `compareRefEnds` (SuffixArrayFuns.cpp:210).
///
/// Decides how two suffixes that hit a spacer compare relative to the
/// `gInsert` insertion boundary (used during sjdb insertion).
#[allow(clippy::if_same_then_else)]
fn compare_ref_ends(map_gen: &Genome, sa_str: u64, g_insert: u64, str_g: bool, str_r: bool) -> i32 {
    if str_g {
        if str_r {
            if sa_str < g_insert {
                1
            } else {
                -1
            }
        } else {
            1
        }
    } else if str_r {
        -1
    } else if g_insert == u64::MAX {
        -1
    } else if sa_str < map_gen.n_genome - g_insert {
        1
    } else {
        -1
    }
}

/// Port of `compareSeqToGenome1` (SuffixArrayFuns.cpp:221).
///
/// Variant of `compareSeqToGenome` that handles `GENOME_spacingChar` in
/// both sequences and `gInsert` boundary for sjdb insertion. Returns
/// `(max_match_length, comp_res)`. Unlike the bool version used by
/// `compareSeqToGenome`, `comp_res` here is a tri-state: negative/zero/positive.
///
/// # Safety
/// `s[0]/s[1]` must extend `≥ S + N` bytes. `map_gen` must have a valid SA
/// and `G` with sufficient padding.
pub unsafe fn compare_seq_to_genome1(
    map_gen: &Genome,
    sa: &PackedArray,
    s2: [&[u8]; 2],
    s: u64,
    n: u64,
    l: u64,
    i_sa: u64,
    dir_r: bool,
    g_insert: u64,
) -> (u64, i32) {
    let sa_str_raw = sa.get(i_sa);
    let dir_g = (sa_str_raw >> map_gen.g_strand_bit) == 0;
    let sa_str = sa_str_raw & map_gen.g_strand_mask;
    let g_base_ptr = unsafe { map_gen.g.as_ptr().add(G_OFFSET) };

    if dir_g {
        let s_ptr = s2[0].as_ptr();
        let s_start = (s + l) as usize;
        let g_start = (sa_str + l) as isize;
        for ii in 0..(n - l) as isize {
            let sv = unsafe { *s_ptr.add(s_start + ii as usize) };
            let gv = unsafe { *g_base_ptr.offset(g_start + ii) };
            if sv != gv {
                return if sv > gv {
                    (ii as u64 + l, 1)
                } else {
                    (ii as u64 + l, -1)
                };
            } else if sv == GENOME_SPACING_CHAR {
                let comp = compare_ref_ends(map_gen, sa_str, g_insert, dir_g, dir_r);
                return (ii as u64 + l, comp);
            }
        }
        (n, 0)
    } else {
        let s_ptr = s2[1].as_ptr();
        let s_start = (s + l) as usize;
        let g_start = (map_gen.n_genome as i64 - 1 - sa_str as i64 - l as i64) as isize;
        for ii in 0..(n - l) as isize {
            let sv = unsafe { *s_ptr.add(s_start + ii as usize) };
            let gv = unsafe { *g_base_ptr.offset(g_start - ii) };
            if sv != gv {
                let mut s1 = sv;
                let mut g1 = gv;
                if s1 < 4 {
                    s1 = 3 - s1;
                }
                if g1 < 4 {
                    g1 = 3 - g1;
                }
                return if s1 > g1 {
                    (ii as u64 + l, 1)
                } else {
                    (ii as u64 + l, -1)
                };
            } else if sv == GENOME_SPACING_CHAR {
                let comp = compare_ref_ends(map_gen, sa_str, g_insert, dir_g, dir_r);
                return (ii as u64 + l, comp);
            }
        }
        (n, 0)
    }
}

/// Port of `suffixArraySearch1` (SuffixArrayFuns.cpp:297).
///
/// Binary search in SA space used by `sjdbBuildIndex`. Returns the first
/// SA index `k` such that `g[SA[k]] > s` (or `SA.len()` if the string is
/// larger than everything in SA). `s[0]/s[1]` are fwd/rev-complement
/// views of the junction sequence.
///
/// # Safety
/// See `compare_seq_to_genome1`.
#[allow(clippy::too_many_arguments)]
pub unsafe fn suffix_array_search1(
    map_gen: &Genome,
    sa: &PackedArray,
    s: [&[u8]; 2],
    ss: u64,
    n: u64,
    g_insert: u64,
    str_r: bool,
    mut i1: u64,
    mut i2: u64,
    mut l: u64,
) -> u64 {
    let (mut l1, comp_res) = compare_seq_to_genome1(map_gen, sa, s, ss, n, l, i1, str_r, g_insert);
    if comp_res < 0 {
        let _ = l1;
        return 0;
    }
    let (mut l2, comp_res) = compare_seq_to_genome1(map_gen, sa, s, ss, n, l, i2, str_r, g_insert);
    if comp_res > 0 {
        let _ = l2;
        return u64::MAX - 1;
    }
    l = std::cmp::min(l1, l2);

    while i1 + 1 < i2 {
        let i3 = median_uint2(i1, i2);
        let (l3, comp_res) = compare_seq_to_genome1(map_gen, sa, s, ss, n, l, i3, str_r, g_insert);
        if l3 == n {
            return i3;
        }
        if comp_res > 0 {
            i1 = i3;
            l1 = l3;
        } else if comp_res < 0 {
            i2 = i3;
            l2 = l3;
        }
        // note: when comp_res == 0 we fall into the else branch above via
        // the `<0` check; preserve C++ semantics by treating equal as
        // "not bigger" → move i2.
        if comp_res == 0 {
            i2 = i3;
            l2 = l3;
        }
        l = std::cmp::min(l1, l2);
    }
    i2
}

/// Port of `funCalcSAi` (SuffixArrayFuns.cpp:397).
///
/// Compute the SAi prefix (length `L`, in 2-bit-per-nt packed form) from a
/// raw genome-numeric byte buffer. Returns `< 0` if the prefix contained a
/// spacer/N before reaching `L`.
pub fn fun_calc_sai(g_seq: &[u8], l: u32) -> i64 {
    let mut ind1: i64 = 0;
    for ii in 0..l as usize {
        let g = g_seq[ii];
        ind1 <<= 2;
        if g > 3 {
            return -1;
        }
        ind1 += g as i64;
    }
    ind1
}

/// Port of `funCompareUintAndSuffixes` (funCompareUintAndSuffixes.cpp:6-40).
///
/// Two-key comparator for `(iSA, iGsj)` tuples used during `qsort` in
/// `sjdbBuildIndex`:
/// 1. Primary: `iSA` ascending.
/// 2. Secondary: suffix of `Gsj` starting at `iGsj`, with `GENOME_spacingChar`
///    acting as an end-of-string marker; identical prefixes up to a
///    spacer fall back to comparing `iGsj`.
pub fn fun_compare_uint_and_suffixes(
    g_sj: &[u8],
    a: &[u64; 2],
    b: &[u64; 2],
) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    match a[0].cmp(&b[0]) {
        Greater => Greater,
        Less => Less,
        Equal => {
            let ga_start = a[1] as usize;
            let gb_start = b[1] as usize;
            let mut ig = 0usize;
            loop {
                let ga = g_sj[ga_start + ig];
                let gb = g_sj[gb_start + ig];
                if ga > gb {
                    return Greater;
                } else if ga < gb {
                    return Less;
                } else if ga == GENOME_SPACING_CHAR {
                    return a[1].cmp(&b[1]);
                }
                ig += 1;
            }
        }
    }
}
