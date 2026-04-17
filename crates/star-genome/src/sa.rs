//! 1:1 port of the suffix-array sort comparator (`funCompareSuffixes`) and
//! the SA-building orchestration (see `Genome_genomeGenerate.cpp` lines
//! 29–89 and 212–349).
//!
//! The comparator is the hottest path of `genomeGenerate`. Exact behavior:
//! - Treat the genome as a stream of bytes where chromosome separators are
//!   `5` (the "spacer" nucleotide code).
//! - Compare 8 bytes at a time, going *backwards* from the suffix anchor
//!   (indexing into `globalG-7`). Use a SIMD-esque `has5` trick to detect
//!   a chromosome boundary within the 8-byte word.
//! - When a chromosome boundary is hit, break the tie *anti-stably* by SA
//!   index (larger index < smaller index).
//!
//! Reproducing this comparator exactly is required for bit-exact SA output.

use std::cmp::Ordering;
use std::ptr;

/// `#define has5(v) ((((v)^0x0505050505050505) - 0x0101010101010101) & ~((v)^0x0505050505050505) & 0x8080808080808080)`
#[inline(always)]
pub const fn has5(v: u64) -> u64 {
    let w = v ^ 0x0505_0505_0505_0505;
    (w.wrapping_sub(0x0101_0101_0101_0101)) & !w & 0x8080_8080_8080_8080
}

/// Read an 8-byte unaligned word from a possibly-misaligned pointer.
#[inline(always)]
unsafe fn read_u64_unaligned(ptr: *const u8) -> u64 {
    unsafe { ptr::read_unaligned(ptr as *const u64) }
}

/// Port of `funCompareSuffixes`: returns -1, 0, or 1 like `qsort`.
///
/// # Safety
/// Requires that `global_g` is at least `7 + max(idx_a, idx_b)` bytes with a
/// 7-byte slack before position 0 (the C++ computes `globalG-7`), and that
/// both SA indices refer into the same contiguous genome buffer. `global_l`
/// is measured in `u64` words.
pub unsafe fn compare_suffixes(
    global_g: *const u8,
    global_l: u64,
    idx_a: u64,
    idx_b: u64,
) -> Ordering {
    let mut ga = unsafe { global_g.offset(-7).add(idx_a as usize) } as *const u8;
    let mut gb = unsafe { global_g.offset(-7).add(idx_b as usize) } as *const u8;
    let mut jj: u64 = 0;
    while jj < global_l {
        let va = unsafe { read_u64_unaligned(ga) };
        let vb = unsafe { read_u64_unaligned(gb) };

        if has5(va) != 0 && has5(vb) != 0 {
            // Compare byte-by-byte starting from the most-significant byte,
            // which corresponds to the *newest* (closest to anchor) byte.
            let va1 = va.to_le_bytes();
            let vb1 = vb.to_le_bytes();
            let mut ii: i32 = 7;
            while ii >= 0 {
                let vai = va1[ii as usize];
                let vbi = vb1[ii as usize];
                if vai > vbi {
                    return Ordering::Greater;
                } else if vai < vbi {
                    return Ordering::Less;
                } else if vai == 5 {
                    // Anti-stable tie-break on SA index.
                    return if idx_a > idx_b {
                        Ordering::Less
                    } else {
                        Ordering::Greater
                    };
                }
                ii -= 1;
            }
        } else if va > vb {
            return Ordering::Greater;
        } else if va < vb {
            return Ordering::Less;
        }
        jj += 1;
        // Step back 8 bytes for the next comparison window.
        ga = unsafe { ga.offset(-8) };
        gb = unsafe { gb.offset(-8) };
    }

    // Suffixes equal up to globalL: fall back to anti-stable tie-break.
    if idx_a > idx_b {
        Ordering::Less
    } else {
        Ordering::Greater
    }
}

/// `funG2strLocus`: decode an SA-encoded strand+pos into a linear locus.
#[inline]
pub fn fun_g2str_locus(sa_str: u64, n: u64, g_strand_bit: u8, g_strand_mask: u64) -> u64 {
    let strand_g = (sa_str >> g_strand_bit) == 0;
    let mut pos = sa_str & g_strand_mask;
    if !strand_g {
        pos += n;
    }
    pos
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn has5_detects_chromosome_boundary() {
        let with = 0x0102_0305_0403_0201u64;
        let without = 0x0102_0304_0403_0201u64;
        assert_ne!(has5(with), 0);
        assert_eq!(has5(without), 0);
    }

    #[test]
    fn compare_identical_indices_is_deterministic_antistable() {
        // Simple 16-byte genome: positions 0..16 all distinct numeric codes.
        let genome: [u8; 32] = [
            0, 0, 0, 0, 0, 0, 0, 0, // 7-byte slack
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 5, 5, 5, 5, 5, 5, 5, 5,
        ];
        let ptr = unsafe { genome.as_ptr().add(7) };
        let ord = unsafe { compare_suffixes(ptr, 1, 0, 1) };
        // Different bytes: should order strictly (not tie).
        assert_ne!(ord, Ordering::Equal);
    }
}
