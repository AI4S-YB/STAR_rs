//! Port of `ReadAlign_stitchPieces.cpp` — orchestrates window creation,
//! per-piece assignment and seed stitching. See the C++ source for inline
//! comments; this file mirrors the control flow 1:1.

use star_core::types::{
    MARKER_NO_GOOD_WINDOW, PC_DIR, PC_I_FRAG, PC_LENGTH, PC_NREP, PC_R_START, PC_SA_END,
    PC_SA_START, UINT_WIN_BIN_MAX, WA_ANCHOR, WC_CHR, WC_G_END, WC_G_START, WC_STR,
};
use star_genome::Genome;
use star_params::parameters::Parameters;

use crate::read_align::ReadAlign;
use crate::stitch_window::{stitch_window_aligns, StitchCtx};
use crate::transcript::Transcript;
use crate::windows::sj_align_split;

impl ReadAlign {
    /// 1:1 port of `ReadAlign::stitchPieces` (ReadAlign_stitchPieces.cpp:12-350).
    ///
    /// # Safety
    /// `r` must have three slices of length ≥ `l_read` for forward / padded /
    /// reverse-complement ASCII bytes. `map_gen.g` must be pre-padded by
    /// `LOAD_L` so negative offsets in `extend_align` are safe.
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn stitch_pieces(
        &mut self,
        p: &Parameters,
        map_gen: &Genome,
        r_buffers: &[&[u8]; 3],
        l_read: u64,
    ) {
        // Zero-out winBin (both strands).
        for strand in 0..2 {
            for cell in self.win_bin[strand].iter_mut() {
                *cell = UINT_WIN_BIN_MAX;
            }
        }

        self.n_w = 0;
        let g_strand_mask = (1u64 << map_gen.g_strand_bit) - 1;

        for i_p in 0..self.n_p as usize {
            if self.pc[i_p][PC_NREP] > p.win_anchor_multimap_nmax {
                continue;
            }
            let a_dir = self.pc[i_p][PC_DIR];
            let a_length = self.pc[i_p][PC_LENGTH];
            for i_sa in self.pc[i_p][PC_SA_START]..=self.pc[i_p][PC_SA_END] {
                let sa_raw = map_gen.sa.get(i_sa);
                let mut a_str = sa_raw >> map_gen.g_strand_bit;
                let mut a1 = sa_raw & g_strand_mask;

                if a_dir == 1 && a_str == 0 {
                    a_str = 1;
                } else if a_dir == 0 && a_str == 1 {
                    a1 = map_gen.n_genome - (a_length + a1);
                } else if a_dir == 1 && a_str == 1 {
                    a_str = 0;
                    a1 = map_gen.n_genome - (a_length + a1);
                }
                if self.revert_strand {
                    a_str = 1 - a_str;
                }

                let mut too_many = false;
                if a1 >= map_gen.sj_g_start {
                    if let Some((a1_d, _a_l_d, a1_a, _a_l_a, _isj)) =
                        sj_align_split(a1, a_length, map_gen)
                    {
                        if self.create_extend_windows_with_align(p, map_gen, a1_d, a_str)
                            == star_core::types::EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS
                        {
                            too_many = true;
                        }
                        if !too_many
                            && self.create_extend_windows_with_align(p, map_gen, a1_a, a_str)
                                == star_core::types::EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS
                        {
                            too_many = true;
                        }
                    }
                } else if self.create_extend_windows_with_align(p, map_gen, a1, a_str)
                    == star_core::types::EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS
                {
                    too_many = true;
                }

                if too_many {
                    break;
                }
            }
        }

        // Flank-extend each window within chromosome limits.
        for i_win in 0..self.n_w as usize {
            if self.wc[i_win][WC_G_START] <= self.wc[i_win][WC_G_END] {
                let mut wb = self.wc[i_win][WC_G_START];
                let chr = self.wc[i_win][WC_CHR];
                let str_ = self.wc[i_win][WC_STR] as usize;
                for _ in 0..p.win_flank_nbins {
                    if wb == 0 {
                        break;
                    }
                    let bin_chr = map_gen.chr_bin[((wb - 1) >> p.win_bin_chr_nbits) as usize];
                    if bin_chr != chr {
                        break;
                    }
                    wb -= 1;
                    if (wb as usize) < self.win_bin[str_].len() {
                        self.win_bin[str_][wb as usize] = i_win as u16;
                    }
                }
                self.wc[i_win][WC_G_START] = wb;

                wb = self.wc[i_win][WC_G_END];
                for _ in 0..p.win_flank_nbins {
                    if wb + 1 >= p.win_bin_n {
                        break;
                    }
                    let bin_chr = map_gen.chr_bin[((wb + 1) >> p.win_bin_chr_nbits) as usize];
                    if bin_chr != chr {
                        break;
                    }
                    wb += 1;
                    if (wb as usize) < self.win_bin[str_].len() {
                        self.win_bin[str_][wb as usize] = i_win as u16;
                    }
                }
                self.wc[i_win][WC_G_END] = wb;
            }
            self.n_wa[i_win] = 0;
            self.wa_lrec[i_win] = 0;
            self.wlast_anchor[i_win] = u64::MAX;
        }
        self.n_wall = self.n_w;

        // Assign pieces to windows.
        for i_p in 0..self.n_p as usize {
            let a_nrep = self.pc[i_p][PC_NREP];
            let a_frag = self.pc[i_p][PC_I_FRAG];
            let a_length = self.pc[i_p][PC_LENGTH];
            let a_dir = self.pc[i_p][PC_DIR];
            let a_anchor = a_nrep <= p.win_anchor_multimap_nmax;

            for ii in 0..self.n_w as usize {
                self.n_wap[ii] = 0;
            }

            for i_sa in self.pc[i_p][PC_SA_START]..=self.pc[i_p][PC_SA_END] {
                let sa_raw = map_gen.sa.get(i_sa);
                let mut a_str = sa_raw >> map_gen.g_strand_bit;
                let mut a1 = sa_raw & g_strand_mask;
                let mut a_rstart = self.pc[i_p][PC_R_START];

                if a_dir == 1 && a_str == 0 {
                    a_str = 1;
                    a_rstart = l_read - (a_length + a_rstart);
                } else if a_dir == 0 && a_str == 1 {
                    a_rstart = l_read - (a_length + a_rstart);
                    a1 = map_gen.n_genome - (a_length + a1);
                } else if a_dir == 1 && a_str == 1 {
                    a_str = 0;
                    a1 = map_gen.n_genome - (a_length + a1);
                }
                if self.revert_strand {
                    a_str = 1 - a_str;
                }

                if a1 >= map_gen.sj_g_start {
                    if let Some((a1_d, a_l_d, a1_a, a_l_a, isj)) =
                        sj_align_split(a1, a_length, map_gen)
                    {
                        self.assign_align_to_window(
                            p, a1_d, a_l_d, a_str, a_nrep, a_frag, a_rstart, a_anchor, isj,
                        );
                        self.assign_align_to_window(
                            p,
                            a1_a,
                            a_l_a,
                            a_str,
                            a_nrep,
                            a_frag,
                            a_rstart + a_l_d,
                            a_anchor,
                            isj,
                        );
                    }
                } else {
                    self.assign_align_to_window(
                        p,
                        a1,
                        a_length,
                        a_str,
                        a_nrep,
                        a_frag,
                        a_rstart,
                        a_anchor,
                        u64::MAX,
                    );
                }
            }
        }

        // Stitch per window.
        if self.tr_all.len() < self.n_w as usize {
            self.tr_all.resize(self.n_w as usize, Vec::new());
        }
        if self.wa_incl.len() < (self.wa.first().map(|v| v.len()).unwrap_or(0)) {
            self.wa_incl
                .resize(self.wa.first().map(|v| v.len()).unwrap_or(32), false);
        }

        let mut tr_best: i32 = 0;
        let mut tr_best_glen: u64 = 0;
        let mut i_w1: u64 = 0;
        for i_w in 0..self.n_w as usize {
            if self.n_wa[i_w] == 0 {
                continue;
            }
            if self.wlast_anchor[i_w] < self.n_wa[i_w] {
                self.wa[i_w][self.wlast_anchor[i_w] as usize][WA_ANCHOR] = 2;
            }
            if self.wa_incl.len() < self.n_wa[i_w] as usize {
                self.wa_incl.resize(self.n_wa[i_w] as usize, false);
            }
            for slot in self.wa_incl.iter_mut().take(self.n_wa[i_w] as usize) {
                *slot = false;
            }

            let mut tr_a = self.tr_init.clone();
            tr_a.chr = self.wc[i_w][WC_CHR];
            tr_a.str_ = self.wc[i_w][WC_STR];
            tr_a.ro_str = if self.revert_strand {
                1 - tr_a.str_
            } else {
                tr_a.str_
            };
            tr_a.max_score = 0;

            // Grow tr_all entry.
            self.tr_all[i_w].clear();
            self.tr_all[i_w].push(tr_a.clone());

            let r_forward = if tr_a.ro_str == 0 {
                r_buffers[0].as_ptr()
            } else {
                r_buffers[2].as_ptr()
            };

            let n_wa_i = self.n_wa[i_w];
            let wa_row = self.wa[i_w].clone();
            let mut w_tr: Vec<Transcript> = vec![tr_a.clone()];
            let mut n_win_tr_i: u64 = 0;

            let read_length = self.read_length;
            let l_read_val = self.l_read;
            let i_read_val = self.i_read;
            let mut max_score_mate = self.max_score_mate;
            let ofmnt = self.out_filter_mismatch_nmax_total;
            let mut ctx = StitchCtx {
                max_score_mate: &mut max_score_mate,
                out_filter_mismatch_nmax_total: ofmnt,
                read_length: &read_length,
                i_read: i_read_val,
            };

            stitch_window_aligns(
                0,
                n_wa_i,
                0,
                &mut self.wa_incl,
                0,
                0,
                tr_a,
                l_read_val,
                &wa_row,
                r_forward,
                map_gen,
                p,
                &mut w_tr,
                &mut n_win_tr_i,
                &mut ctx,
            );
            self.max_score_mate = max_score_mate;

            if n_win_tr_i == 0 {
                continue;
            }
            self.n_win_tr[i_w1 as usize] = n_win_tr_i;
            self.tr_all[i_w1 as usize] = w_tr;
            let best = &self.tr_all[i_w1 as usize][0];
            if best.max_score > tr_best
                || (best.max_score == tr_best && best.g_length < tr_best_glen)
            {
                tr_best = best.max_score;
                tr_best_glen = best.g_length;
                self.tr_best = i_w1;
            }
            i_w1 += 1;
        }

        self.n_w = i_w1;

        if self.tr_best == u64::MAX || tr_best == 0 {
            self.map_marker = MARKER_NO_GOOD_WINDOW;
            self.n_w = 0;
        }
    }
}
