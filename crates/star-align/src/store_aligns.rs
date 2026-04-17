//! Port of `ReadAlign_storeAligns.cpp`.
//!
//! Keeps the `OPTIM_STOREaligns_SIMPLE` fast-path behavior of the upstream
//! code.

use star_core::types::{
    EXIT_CODE_RUNTIME, PC_DIR, PC_I_FRAG, PC_LENGTH, PC_NREP, PC_R_START, PC_SA_END, PC_SA_START,
    PC_SIZE,
};
use star_params::parameters::Parameters;

use crate::read_align::{PieceCoord, ReadAlign};

impl ReadAlign {
    /// 1:1 port of `ReadAlign::storeAligns` (`ReadAlign_storeAligns.cpp:10-160`).
    pub fn store_aligns(
        &mut self,
        p: &Parameters,
        i_dir: u64,
        shift: u64,
        nrep: u64,
        mut l: u64,
        ind_start_end: [u64; 2],
        i_frag: u64,
    ) -> anyhow::Result<()> {
        if nrep > p.seed_multimap_nmax {
            if nrep < self.mult_nmin || self.mult_nmin == 0 {
                self.mult_nmin = nrep;
                self.mult_nmin_l = l;
            }
            return Ok(());
        }

        let idx = if nrep == 1 { 0 } else { 1 };
        self.n_um[idx] += nrep;
        self.n_a += nrep;

        let r_start = if i_dir == 0 { shift } else { shift + 1 - l };

        // Ensure PC has capacity for one extra slot (insertion point may be at end).
        if self.pc.len() <= self.n_p as usize + 1 {
            let grow_to = (self.n_p as usize + 2).max(32);
            self.pc.resize(grow_to, [0u64; PC_SIZE]);
        }

        let mut i_p: i64 = self.n_p as i64 - 1;
        while i_p >= 0 {
            if self.pc[i_p as usize][0] <= r_start {
                if self.pc[i_p as usize][PC_R_START] == r_start
                    && self.pc[i_p as usize][PC_LENGTH] < l
                {
                    i_p -= 1;
                    continue;
                }
                if self.pc[i_p as usize][PC_R_START] == r_start
                    && self.pc[i_p as usize][PC_LENGTH] == l
                {
                    return Ok(());
                }
                break;
            }
            i_p -= 1;
        }
        let i_p = (i_p + 1) as usize;

        for ii in (i_p..self.n_p as usize).rev() {
            let prev: PieceCoord = self.pc[ii];
            self.pc[ii + 1] = prev;
        }

        self.n_p += 1;
        if self.n_p > p.seed_per_read_nmax {
            anyhow::bail!(
                "EXITING because of FATAL error (exit code {}): too many pieces per read — raise --seedPerReadNmax",
                EXIT_CODE_RUNTIME
            );
        }

        self.pc[i_p][PC_R_START] = r_start;
        self.pc[i_p][PC_LENGTH] = l;
        self.pc[i_p][PC_DIR] = i_dir;
        self.pc[i_p][PC_NREP] = nrep;
        self.pc[i_p][PC_SA_START] = ind_start_end[0];
        self.pc[i_p][PC_SA_END] = ind_start_end[1];
        self.pc[i_p][PC_I_FRAG] = i_frag;

        if l < self.stored_lmin {
            l = self.stored_lmin;
        }

        if nrep == 1 {
            if l > self.uniq_lmax {
                self.uniq_lmax = l;
                self.uniq_lmax_ind = self.n_p - 1;
            }
        } else {
            if nrep < self.mult_nmin || self.mult_nmin == 0 {
                self.mult_nmin = nrep;
                self.mult_nmin_l = l;
            }
            if l > self.mult_lmax {
                self.mult_lmax = l;
                self.mult_lmax_n = nrep;
            }
            if nrep > self.mult_nmax {
                self.mult_nmax = nrep;
                self.mult_nmax_l = l;
            }
        }
        Ok(())
    }
}
