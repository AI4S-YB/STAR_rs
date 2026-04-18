//! 1:1 port of `ReadAlign::mapOneRead` (ReadAlign_mapOneRead.cpp).
//!
//! Drives the seed-finding passes then, if enough seeds are found, calls
//! `stitch_pieces` to build the final set of candidate transcripts.

use star_core::seq::quality_split;
use star_core::types::{
    MARKER_ALL_PIECES_EXCEED_SEED_MULTIMAP_NMAX, MARKER_NO_GOOD_PIECES, MARKER_READ_TOO_SHORT,
};
use star_genome::Genome;
use star_params::parameters::Parameters;

use crate::read_align::ReadAlign;
use crate::seed::max_mappable_length_2strands;

impl ReadAlign {
    /// 1:1 port of `mapOneRead()`.
    ///
    /// `r_buffers` must hold 3 ASCII-numeric-encoded buffers:
    ///   `[0]` forward, `[1]` complement (non-reversed), `[2]` reverse-complement.
    /// All three must be at least `l_read` bytes long.
    ///
    /// # Safety
    /// The borrow split between `self` and the read buffers is enforced by the
    /// caller (no self-references into the buffers). `map_gen` must be fully
    /// loaded including `LOAD_L` padding.
    pub unsafe fn map_one_read(
        &mut self,
        p: &Parameters,
        map_gen: &Genome,
        r_buffers: &[&[u8]; 3],
    ) -> i32 {
        // The 2nd read is always on the opposite strand — it is reversed
        // outside this function, so revertStrand stays false.
        self.revert_strand = false;

        if self.l_read > 0 {
            self.ensure_split_capacity(p.max_n_split as usize);
            self.n_split = quality_split(
                &r_buffers[0][..self.l_read as usize],
                p.max_n_split,
                p.seed_split_min,
                &mut self.split_r,
            );
        } else {
            self.n_split = 0;
        }

        self.reset_n();

        // Prime the template transcript (used as a starting template for each
        // window).
        self.tr_init.reset();
        self.tr_init.chr = 0;
        self.tr_init.str_ = 0;
        self.tr_init.ro_str = 0;
        self.tr_init.c_start = 0;
        self.tr_init.g_length = 0;
        self.tr_init.i_read = self.i_read;
        self.tr_init.l_read = self.l_read;
        self.tr_init.n_exons = 0;
        self.tr_init.read_length_pair_original = self.read_length_pair_original;
        self.tr_init.read_n_mates = self.read_n_mates as u64;
        self.tr_init.read_length_pair_original = self.read_length_pair_original;

        self.tr_best = u64::MAX;

        let seed_search_start_lmax = std::cmp::min(
            p.seed_search_start_lmax,
            (p.seed_search_start_lmax_over_lread * (self.l_read.saturating_sub(1)) as f64) as u64,
        );

        for i_p in 0..self.n_split as usize {
            let piece_len = self.split_r[i_p][1];
            let piece_start = self.split_r[i_p][0];
            let i_frag = self.split_r[i_p][2] as u32;

            let n_start = if p.seed_search_start_lmax > 0 && seed_search_start_lmax < piece_len {
                piece_len / seed_search_start_lmax + 1
            } else {
                1
            };
            let l_start = piece_len / n_start;
            let mut flag_dir_map = true;

            for i_dir in 0..2u32 {
                for istart in 0..n_start {
                    if !(flag_dir_map || istart > 0) {
                        continue;
                    }
                    let mut l_mapped: u64 = 0;
                    while istart * l_start + l_mapped + p.seed_map_min < piece_len {
                        let shift = if i_dir == 0 {
                            piece_start + istart * l_start + l_mapped
                        } else {
                            piece_start + piece_len - istart * l_start - 1 - l_mapped
                        };
                        let seed_length = piece_len - l_mapped - istart * l_start;
                        let (_nrep, l, seeds) = max_mappable_length_2strands(
                            map_gen,
                            p,
                            &[r_buffers[0], r_buffers[1]],
                            shift,
                            seed_length,
                            i_dir,
                            0,
                            map_gen.n_sa - 1,
                            i_frag,
                        );
                        for seed in &seeds {
                            self.store_aligns(
                                p,
                                seed.i_dir as u64,
                                seed.shift,
                                seed.nrep,
                                seed.max_l,
                                seed.ind_start_end,
                                seed.i_frag as u64,
                            )
                            .ok();
                        }
                        if i_dir == 0 && istart == 0 && l_mapped == 0 && shift + l == piece_len {
                            flag_dir_map = false;
                        }
                        // C++ does `Lmapped += L;` unconditionally — if L==0 the
                        // original code infinite-loops. That edge case never
                        // triggers on real reads (every position in a real
                        // genome matches something in the SA), but our small
                        // test fixtures can hit it. Break instead of hanging;
                        // no stored seeds are lost because L==0 implies no hit
                        // was returned for this shift.
                        if l == 0 {
                            break;
                        }
                        l_mapped += l;
                    }
                }

                if p.seed_search_lmax > 0 {
                    for istart in 0..n_start {
                        let shift = if i_dir == 0 {
                            piece_start + istart * l_start
                        } else {
                            piece_start + piece_len - istart * l_start - 1
                        };
                        let seed_length = if i_dir == 0 {
                            std::cmp::min(p.seed_search_lmax, piece_start + piece_len - shift)
                        } else {
                            std::cmp::min(p.seed_search_lmax, shift + 1)
                        };
                        let (_nrep, _l, seeds) = max_mappable_length_2strands(
                            map_gen,
                            p,
                            &[r_buffers[0], r_buffers[1]],
                            shift,
                            seed_length,
                            i_dir,
                            0,
                            map_gen.n_sa - 1,
                            i_frag,
                        );
                        for seed in &seeds {
                            self.store_aligns(
                                p,
                                seed.i_dir as u64,
                                seed.shift,
                                seed.nrep,
                                seed.max_l,
                                seed.ind_start_end,
                                seed.i_frag as u64,
                            )
                            .ok();
                        }
                    }
                }
            }
        }

        if (self.l_read as i32) < p.out_filter_match_nmin {
            self.map_marker = MARKER_READ_TOO_SHORT;
            self.tr_init.r_length = 0;
            self.n_w = 0;
        } else if self.n_split == 0 {
            self.map_marker = MARKER_NO_GOOD_PIECES;
            self.tr_init.r_length = self.split_r.first().map(|r| r[1]).unwrap_or(0);
            self.n_w = 0;
        } else if self.n_split > 0 && self.n_a == 0 {
            self.map_marker = MARKER_ALL_PIECES_EXCEED_SEED_MULTIMAP_NMAX;
            self.tr_init.r_length = self.mult_nmin_l;
            self.n_w = 0;
        } else if self.n_split > 0 && self.n_a > 0 {
            self.stitch_pieces(p, map_gen, r_buffers, self.l_read);
        }

        0
    }

    fn ensure_split_capacity(&mut self, n: usize) {
        if self.split_r.len() < n.max(1) {
            self.split_r.resize(n.max(1), [0u64; 3]);
        }
    }
}
