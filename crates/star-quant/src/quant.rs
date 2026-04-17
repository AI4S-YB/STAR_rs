//! 1:1 port of `Quantifications.{h,cpp}` + the gene-count portion of
//! `Transcriptome_geneCountsAddAlign.cpp`.
//!
//! Per-chunk counter arrays + the `geneCountsAddAlign` classifier. Merge
//! via [`Quantifications::add_quants`].

use star_align::Transcript;
use star_core::{EX_G, EX_L};

use crate::transcriptome::ExG;

/// Gene counts for the three "type" slots (unstranded / forward / reverse).
///
/// Port of the nested `geneCounts` struct in `Quantifications.h`.
#[derive(Debug, Clone)]
pub struct GeneCounts {
    /// Always 3 (unstranded, same-strand, reverse-strand).
    pub n_type: usize,
    pub n_ge: u32,
    pub c_multi: u64,
    /// Per-type count of reads that could not be assigned unambiguously
    /// to a single gene.
    pub c_ambig: [u64; 3],
    /// Per-type count of reads that did not overlap any exon.
    pub c_none: [u64; 3],
    /// `g_count[itype][ig]` — per-gene, per-type counts. Layout mirrors
    /// C++ `uintQ **gCount`; we flatten to `Vec<Vec<u64>>` for
    /// simplicity.
    pub g_count: [Vec<u64>; 3],
}

impl GeneCounts {
    pub fn new(n_ge: u32) -> Self {
        Self {
            n_type: 3,
            n_ge,
            c_multi: 0,
            c_ambig: [0; 3],
            c_none: [0; 3],
            g_count: [
                vec![0; n_ge as usize],
                vec![0; n_ge as usize],
                vec![0; n_ge as usize],
            ],
        }
    }
}

/// Port of `Quantifications` (minus the TranscriptomeSAM parts, which are
/// not in M7.3 scope).
#[derive(Debug, Clone)]
pub struct Quantifications {
    pub gene_counts: GeneCounts,
}

impl Quantifications {
    pub fn new(n_ge: u32) -> Self {
        Self {
            gene_counts: GeneCounts::new(n_ge),
        }
    }

    /// Port of `Quantifications::addQuants` (lines 25-37 of Quantifications.cpp).
    pub fn add_quants(&mut self, other: &Quantifications) {
        self.gene_counts.c_multi += other.gene_counts.c_multi;
        for it in 0..self.gene_counts.n_type {
            self.gene_counts.c_ambig[it] += other.gene_counts.c_ambig[it];
            self.gene_counts.c_none[it] += other.gene_counts.c_none[it];
            for ig in 0..self.gene_counts.n_ge as usize {
                self.gene_counts.g_count[it][ig] += other.gene_counts.g_count[it][ig];
            }
        }
    }
}

/// Port of `Transcriptome::geneCountsAddAlign`
/// (`Transcriptome_geneCountsAddAlign.cpp` lines 4-60).
///
/// Takes the set of best alignments (`aligns[0..n_a]`), classifies the
/// read for each of the three strand slots (unstranded / same / reverse),
/// and increments the appropriate counters on `quants`.
///
/// `gene1` is filled with the selected gene id per slot (`-1` = no feature,
/// `-2` = ambiguous); it mirrors the `vector<int32> &gene1` out-parameter
/// in C++ and is used by the Solo code path. For the pure GeneCounts case,
/// callers can pass a scratch `[-1; 3]` array.
pub fn gene_counts_add_align(
    quants: &mut Quantifications,
    ex_g: &ExG,
    aligns: &[&Transcript],
    gene1: &mut [i32; 3],
) {
    gene1.fill(-1);
    let n_a = aligns.len();

    if n_a > 1 {
        quants.gene_counts.c_multi += 1;
        return;
    }
    if n_a == 0 {
        return;
    }

    let a = aligns[0];

    // Scan all blocks (exons) of the alignment from last to first.
    for ib in (0..a.n_exons as usize).rev() {
        let g_end = a.exons[ib][EX_G] + a.exons[ib][EX_L] - 1; // end of block

        // Binary search: find the largest i with ex_g.s[i] <= g_end. Port of
        // `binarySearch1a<uint64>` from STAR/source/serviceFuns.cpp.
        let mut e1: i64 = binary_search1a_u64(g_end, &ex_g.s);

        while e1 >= 0 && ex_g.e_max[e1 as usize] >= a.exons[ib][EX_G] {
            if ex_g.e[e1 as usize] >= a.exons[ib][EX_G] {
                let str1 = ex_g.str_[e1 as usize] as u32 - 1; // 0 = +, 1 = -, 2 = unstranded
                for itype in 0..3usize {
                    // Strand filter (match C++: `if (itype==1 && a.Str!=str1 && str1<2) continue;`)
                    if itype == 1 && (a.str_ as u32) != str1 && str1 < 2 {
                        continue;
                    }
                    if itype == 2 && (a.str_ as u32) == str1 && str1 < 2 {
                        continue;
                    }
                    match gene1[itype] {
                        -1 => gene1[itype] = ex_g.g[e1 as usize] as i32,
                        -2 => continue,
                        cur if cur != ex_g.g[e1 as usize] as i32 => gene1[itype] = -2,
                        _ => {}
                    }
                }
            }
            e1 -= 1;
        }
    }

    for itype in 0..3usize {
        match gene1[itype] {
            -1 => quants.gene_counts.c_none[itype] += 1,
            -2 => quants.gene_counts.c_ambig[itype] += 1,
            g => quants.gene_counts.g_count[itype][g as usize] += 1,
        }
    }
}

/// Port of `binarySearch1a<uint64>` — finds the largest index `i` such
/// that `arr[i] <= key`. Returns `-1` if all `arr[i] > key`.
///
/// C++ source: `STAR/source/serviceFuns.cpp`. The array is assumed to be
/// sorted non-decreasing; if multiple elements are equal to `key`, the
/// index of the last is returned.
fn binary_search1a_u64(key: u64, arr: &[u64]) -> i64 {
    if arr.is_empty() || arr[0] > key {
        return -1;
    }
    let mut lo: i64 = 0;
    let mut hi: i64 = arr.len() as i64 - 1;
    // Invariant: arr[lo] <= key < arr[hi+1] (sentinel past-the-end).
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if arr[mid as usize] > key {
            hi = mid - 1;
        } else {
            lo = mid;
        }
    }
    // Final refinement
    if hi > lo && arr[hi as usize] <= key {
        hi
    } else {
        lo
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binsearch_basic() {
        let a = [10u64, 20, 30, 40];
        assert_eq!(binary_search1a_u64(5, &a), -1);
        assert_eq!(binary_search1a_u64(10, &a), 0);
        assert_eq!(binary_search1a_u64(15, &a), 0);
        assert_eq!(binary_search1a_u64(20, &a), 1);
        assert_eq!(binary_search1a_u64(39, &a), 2);
        assert_eq!(binary_search1a_u64(40, &a), 3);
        assert_eq!(binary_search1a_u64(1000, &a), 3);
    }

    #[test]
    fn binsearch_duplicates() {
        let a = [5u64, 5, 5, 10];
        assert_eq!(binary_search1a_u64(5, &a), 2);
    }
}
