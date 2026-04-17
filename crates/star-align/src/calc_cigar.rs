//! 1:1 port of `ReadAlign::calcCIGAR` (ReadAlign_calcCIGAR.cpp).

use std::fmt::Write as _;

use star_core::types::{EX_G, EX_I_FRAG, EX_L, EX_R};

use crate::read_align::ReadAlign;
use crate::transcript::Transcript;

impl ReadAlign {
    /// 1:1 port of `calcCIGAR(trOut, nMates, iExMate, leftMate)`.
    /// Fills `self.mates_cigar` with one CIGAR string per mate.
    pub fn calc_cigar(
        &mut self,
        tr_out: &Transcript,
        n_mates: u64,
        i_ex_mate: u64,
        left_mate: u64,
    ) {
        self.mates_cigar.clear();
        for imate in 0..n_mates as usize {
            let i_ex1 = if imate == 0 { 0u64 } else { i_ex_mate + 1 };
            let i_ex2 = if imate == 0 {
                i_ex_mate
            } else {
                tr_out.n_exons - 1
            };
            let mate = tr_out.exons[i_ex1 as usize][EX_I_FRAG];
            let str_ = tr_out.str_;

            let mut cigar = String::new();

            let trim_l = if str_ == 0 && mate == 0 {
                self.clip_mates[mate as usize][0].clipped_n
            } else if str_ == 0 && mate == 1 {
                self.clip_mates[mate as usize][1].clipped_n
            } else if str_ == 1 && mate == 0 {
                self.clip_mates[mate as usize][1].clipped_n
            } else {
                self.clip_mates[mate as usize][0].clipped_n
            };

            let ex_r = tr_out.exons[i_ex1 as usize][EX_R];
            let sub = if ex_r < self.read_length[left_mate as usize] {
                0
            } else {
                self.read_length[left_mate as usize] + 1
            };
            let trim_l1 = trim_l + ex_r - sub;
            if trim_l1 > 0 {
                write!(&mut cigar, "{}S", trim_l1).unwrap();
            }

            for ii in i_ex1..=i_ex2 {
                let ii_u = ii as usize;
                if ii > i_ex1 {
                    let prev = &tr_out.exons[ii_u - 1];
                    let cur = &tr_out.exons[ii_u];
                    let gap_g = cur[EX_G] - (prev[EX_G] + prev[EX_L]);
                    let gap_r = cur[EX_R] - prev[EX_R] - prev[EX_L];
                    if gap_r > 0 {
                        write!(&mut cigar, "{}I", gap_r).unwrap();
                    }
                    if tr_out.canon_sj[ii_u - 1] >= 0 || tr_out.sj_annot[ii_u - 1] == 1 {
                        write!(&mut cigar, "{}N", gap_g).unwrap();
                    } else if gap_g > 0 {
                        write!(&mut cigar, "{}D", gap_g).unwrap();
                    }
                }
                write!(&mut cigar, "{}M", tr_out.exons[ii_u][EX_L]).unwrap();
            }

            let lr_left = self.read_length[left_mate as usize];
            let lr_orig_left = self.read_length_original[left_mate as usize];
            let lr_orig_mate = self.read_length_original[mate as usize];
            let ex1_r = tr_out.exons[i_ex1 as usize][EX_R];
            let ex2_r = tr_out.exons[i_ex2 as usize][EX_R];
            let ex2_l = tr_out.exons[i_ex2 as usize][EX_L];
            let base = if ex1_r < lr_left {
                lr_orig_left
            } else {
                lr_left + 1 + lr_orig_mate
            };
            let trim_r1 = base.saturating_sub(ex2_r + ex2_l + trim_l);
            if trim_r1 > 0 {
                write!(&mut cigar, "{}S", trim_r1).unwrap();
            }

            self.mates_cigar.push(cigar);
        }
    }
}
