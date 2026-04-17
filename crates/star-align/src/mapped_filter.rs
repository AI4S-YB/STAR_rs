//! 1:1 port of `ReadAlign::mappedFilter` (ReadAlign_mappedFilter.cpp).

use star_params::parameters::Parameters;
use star_stats::stats::Stats;

use crate::read_align::ReadAlign;

impl ReadAlign {
    /// 1:1 port of `mappedFilter()`. Returns `true` if the read passes
    /// (unmapType = -1); otherwise updates `stats_ra` counters and `unmap_type`.
    pub fn mapped_filter(&mut self, p: &Parameters, stats_ra: &mut Stats) -> bool {
        self.unmap_type = -1;
        let best = if self.n_w > 0 {
            self.tr_all
                .get(self.tr_best as usize)
                .and_then(|v| v.first())
        } else {
            None
        };

        if self.n_w == 0 || best.is_none() {
            stats_ra.unmapped_other += 1;
            self.unmap_type = 0;
            return false;
        }
        let best = best.unwrap();
        let l_read = self.l_read;

        // C++ `mappedFilter` casts the float-scaled threshold to int
        // (`(intScore)(...)` / `(uint)(...)` in ReadAlign_mappedFilter.cpp:8-9)
        // before comparing. Preserve that truncation so that a score equal
        // to `floor(frac * (Lread-1))` passes, e.g. score 199 vs 0.66*302=199.32.
        let l_read_m1 = l_read.saturating_sub(1) as f64;
        let score_min_frac = (p.out_filter_score_min_over_lread * l_read_m1) as i32;
        let match_min_frac = (p.out_filter_match_nmin_over_lread * l_read_m1) as u64;
        if best.max_score < p.out_filter_score_min
            || best.max_score < score_min_frac
            || (best.n_match as i32) < p.out_filter_match_nmin
            || best.n_match < match_min_frac
        {
            stats_ra.unmapped_short += 1;
            self.unmap_type = 1;
            return false;
        }

        if best.n_mm > self.out_filter_mismatch_nmax_total
            || (best.n_mm as f64 / best.r_length as f64) > p.out_filter_mismatch_nover_lmax
        {
            stats_ra.unmapped_mismatch += 1;
            self.unmap_type = 2;
            return false;
        }

        if self.n_tr > p.out_filter_multimap_nmax as u64 {
            stats_ra.unmapped_multi += 1;
            self.unmap_type = 3;
            return false;
        }

        true
    }
}
