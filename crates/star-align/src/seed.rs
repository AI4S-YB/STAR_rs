//! 1:1 port of `ReadAlign::maxMappableLength2strands`
//! (ReadAlign_maxMappableLength2strands.cpp).
//!
//! The C++ method hangs off of `ReadAlign` and mutates many instance fields.
//! We split that into a pure seed-search function here that returns the
//! [`SeedResult`] list; the caller ([`ReadAlignChunk`] / `ReadAlign::mapOneRead`)
//! decides how to integrate the results into its piece table (`PC` array).

use star_core::UInt;
use star_genome::genome::Genome;
use star_genome::sa_funs::{compare_seq_to_genome, max_mappable_length};
use star_params::parameters::Parameters;

/// One seed hit produced by [`max_mappable_length_2strands`]. Corresponds to
/// a single `storeAligns` call in the C++ code.
#[derive(Debug, Clone, Copy)]
pub struct SeedHit {
    pub i_dir: u32,
    pub shift: u64,
    pub nrep: u64,
    pub max_l: u64,
    pub ind_start_end: [u64; 2],
    pub i_frag: u32,
}

/// Port of `maxMappableLength2strands`.
///
/// Returns the tuple `(Nrep, max_l_best, seeds)` where `seeds` is the list
/// of hits that the C++ code would pass to `storeAligns`.
///
/// # Safety
/// The caller must guarantee that `read[0]` / `read[1]` are valid for reads
/// covering `piece_start_in .. piece_start_in + piece_length_in` (forward) or
/// the corresponding reverse range. `genome` must be fully loaded.
pub unsafe fn max_mappable_length_2strands(
    map_gen: &Genome,
    p: &Parameters,
    read: &[&[u8]; 2],
    piece_start_in: u64,
    piece_length_in: u64,
    i_dir: u32,
    mut i_sa1: u64,
    _i_sa2: u64,
    i_frag: u32,
) -> (u64, u64, Vec<SeedHit>) {
    let sparse_d = p.p_ge.g_sa_sparse_d.max(1);
    let dir_r = i_dir == 0;
    let index_nbases = p.p_ge.g_sa_index_nbases;

    let mut nrep_final: u64 = 0;
    let mut max_l_best: u64 = 0;
    let mut max_l_all = vec![0u64; sparse_d as usize];
    let mut nrep_all = vec![0u64; sparse_d as usize];
    let mut ind_start_end_all = vec![[0u64; 2]; sparse_d as usize];

    let limit = std::cmp::min(piece_length_in, sparse_d);

    for i_dist in 0..limit {
        let piece_length = piece_length_in - i_dist;
        let l_max = std::cmp::min(index_nbases, piece_length);
        let mut ind1: u64 = 0;
        let piece_start: u64 = if dir_r {
            let start = piece_start_in + i_dist;
            for ii in 0..l_max {
                ind1 <<= 2;
                ind1 += read[0][(start + ii) as usize] as u64;
            }
            start
        } else {
            let start = piece_start_in - i_dist;
            for ii in 0..l_max {
                ind1 <<= 2;
                ind1 += 3 - read[0][(start - ii) as usize] as u64;
            }
            start
        };

        let mut l_ind = l_max;
        while l_ind > 0 {
            let sa_idx = map_gen.genome_sa_index_start[(l_ind - 1) as usize] + ind1;
            i_sa1 = map_gen.sai.get(sa_idx);
            if (i_sa1 & map_gen.sai_mark_absent_mask_c) == 0 {
                break;
            }
            l_ind -= 1;
            ind1 >>= 2;
        }

        let mut i_sa2_good = true;
        let i_sa2;
        let base = map_gen.genome_sa_index_start[(l_ind - 1) as usize] + ind1 + 1;
        if base < map_gen.genome_sa_index_start[l_ind as usize] {
            let v = map_gen.sai.get(base);
            if (v & map_gen.sai_mark_absent_mask_c) == 0 {
                i_sa2 = (v & map_gen.sai_mark_n_mask) - 1;
            } else {
                i_sa2 = map_gen.n_sa - 1;
                i_sa2_good = false;
            }
        } else {
            i_sa2 = map_gen.n_sa - 1;
            i_sa2_good = false;
        }

        let i_sa1_no_n = (i_sa1 & map_gen.sai_mark_n_mask_c) == 0;

        let (nrep, max_l, ind_start_end);
        if l_ind < index_nbases && i_sa1_no_n && i_sa2_good {
            let a = i_sa1;
            let b = i_sa2;
            nrep = b - a + 1;
            max_l = l_ind;
            ind_start_end = [a, b];
        } else if i_sa1 == i_sa2 && i_sa1_no_n && i_sa2_good {
            nrep = 1;
            let (ml, _) = compare_seq_to_genome(
                map_gen,
                read,
                piece_start,
                piece_length,
                l_ind,
                i_sa1,
                dir_r,
            );
            max_l = ml;
            ind_start_end = [i_sa1, i_sa2];
        } else {
            let start_l = if i_sa2_good && i_sa1_no_n { l_ind } else { 0 };
            let (n, l, ind) = max_mappable_length(
                map_gen,
                read,
                piece_start,
                piece_length,
                i_sa1 & map_gen.sai_mark_n_mask,
                i_sa2,
                dir_r,
                start_l,
            );
            nrep = n;
            max_l = l;
            ind_start_end = ind;
        }

        if max_l + i_dist > max_l_best {
            max_l_best = max_l + i_dist;
        }
        nrep_all[i_dist as usize] = nrep;
        max_l_all[i_dist as usize] = max_l;
        ind_start_end_all[i_dist as usize] = ind_start_end;
        nrep_final = nrep;
    }

    let mut seeds = Vec::new();
    for i_dist in 0..limit {
        if max_l_all[i_dist as usize] + i_dist == max_l_best {
            let shift = if dir_r {
                piece_start_in + i_dist
            } else {
                piece_start_in - i_dist
            };
            seeds.push(SeedHit {
                i_dir,
                shift,
                nrep: nrep_all[i_dist as usize],
                max_l: max_l_all[i_dist as usize],
                ind_start_end: ind_start_end_all[i_dist as usize],
                i_frag,
            });
        }
    }

    let _ = nrep_final;
    let nrep_final = seeds.last().map(|s| s.nrep).unwrap_or(0);
    (nrep_final, max_l_best, seeds)
}

#[allow(dead_code)]
pub(crate) fn _touch_uint(_: UInt) {}
