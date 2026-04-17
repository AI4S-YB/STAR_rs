//! 1:1 port of `ReadAlign::multMapSelect` (ReadAlign_multMapSelect.cpp).

use star_core::types::{EXIT_CODE_BUG, EXIT_CODE_PARAMETER, MAX_N_MULTMAP};
use star_genome::Genome;
use star_params::parameters::Parameters;

use crate::read_align::ReadAlign;
use crate::rng::rng_uniform_real_0_to_1;

impl ReadAlign {
    /// 1:1 port of `multMapSelect()`.
    pub fn mult_map_select(&mut self, p: &Parameters, map_gen: &Genome) -> anyhow::Result<()> {
        self.n_tr = 0;
        self.tr_mult_array.clear();
        if self.n_w == 0 {
            return Ok(());
        }

        self.max_score = -10 * self.l_read as i32;
        for i_w in 0..self.n_w as usize {
            if let Some(best) = self.tr_all.get(i_w).and_then(|v| v.first()) {
                if self.max_score < best.max_score {
                    self.max_score = best.max_score;
                }
            }
        }

        let tr_best_score = self
            .tr_all
            .get(self.tr_best as usize)
            .and_then(|v| v.first())
            .map(|t| t.max_score)
            .unwrap_or(i32::MIN);
        if self.max_score != tr_best_score {
            anyhow::bail!(
                "EXITING because of FATAL ERROR (exit code {}): BUG: maxScore!=trBest->maxScore in multMapSelect",
                EXIT_CODE_BUG
            );
        }

        // Gather qualifying transcripts.
        for i_w in 0..self.n_w as usize {
            let n_win_tr = self.n_win_tr.get(i_w).copied().unwrap_or(0) as usize;
            for i_tr in 0..n_win_tr {
                let tr = self.tr_all[i_w][i_tr].clone();
                if (tr.max_score + p.out_filter_multimap_score_range as i32) >= self.max_score {
                    if self.n_tr as usize == MAX_N_MULTMAP {
                        anyhow::bail!(
                            "EXITING because of Fatal ERROR (exit code {}): number of alignments exceeds MAX_N_MULTMAP",
                            EXIT_CODE_PARAMETER
                        );
                    }
                    let mut chosen = tr;
                    let first = &self.tr_all[i_w][0];
                    chosen.chr = first.chr;
                    chosen.str_ = first.str_;
                    chosen.ro_str = first.ro_str;
                    self.tr_mult_array.push(chosen);
                    self.n_tr += 1;
                }
            }
        }

        if self.n_tr > p.out_filter_multimap_nmax as u64 || self.n_tr == 0 {
            return Ok(());
        }

        for tr in self.tr_mult_array.iter_mut() {
            tr.ro_start = if tr.ro_str == 0 {
                tr.r_start
            } else {
                self.l_read - tr.r_start - tr.r_length
            };
            tr.c_start = tr.g_start - map_gen.chr_start[tr.chr as usize];
        }

        if self.n_tr == 1 {
            self.tr_mult_array[0].primary_flag = true;
        } else {
            let mut nbest: usize = 0;
            let u_sentinel = u64::MAX;
            let out_sam_mult_nmax_set = (p.out_sam_multi_nmax as i64) != -1;
            if p.out_multimapper_order_random || out_sam_mult_nmax_set {
                let n_tr = self.n_tr as usize;
                let max_score = self.max_score;
                for itr in 0..n_tr {
                    if self.tr_mult_array[itr].max_score == max_score {
                        self.tr_mult_array.swap(itr, nbest);
                        nbest += 1;
                    }
                }
            }
            if p.out_multimapper_order_random {
                if nbest >= 2 {
                    for itr in (1..nbest).rev() {
                        let r = rng_uniform_real_0_to_1(&mut self.rng_mult_order);
                        let rand1 = (r * itr as f64 + 0.5) as usize;
                        self.tr_mult_array.swap(itr, rand1);
                    }
                }
                let n_tr = self.n_tr as usize;
                if n_tr.saturating_sub(nbest) >= 2 {
                    for itr in (1..(n_tr - nbest)).rev() {
                        let r = rng_uniform_real_0_to_1(&mut self.rng_mult_order);
                        let rand1 = (r * itr as f64 + 0.5) as usize;
                        self.tr_mult_array.swap(nbest + itr, nbest + rand1);
                    }
                }
            }
            if p.out_sam_primary_flag == "AllBestScore" {
                let max_score = self.max_score;
                for tr in self.tr_mult_array.iter_mut() {
                    if tr.max_score == max_score {
                        tr.primary_flag = true;
                    }
                }
            } else if p.out_multimapper_order_random || out_sam_mult_nmax_set {
                self.tr_mult_array[0].primary_flag = true;
            } else {
                let _ = u_sentinel;
                let best_idx = self
                    .tr_mult_array
                    .iter()
                    .position(|t| t.max_score == self.max_score)
                    .unwrap_or(0);
                self.tr_mult_array[best_idx].primary_flag = true;
            }
        }

        Ok(())
    }
}
