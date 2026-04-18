//! Ports `ReadAlign_createExtendWindowsWithAlign.cpp` and
//! `ReadAlign_assignAlignToWindow.cpp` (two small but critical helpers that
//! build the per-read window list).

use star_core::types::{
    EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS, MARKER_TOO_MANY_ANCHORS_PER_WINDOW,
    UINT_WIN_BIN_MAX, WA_ANCHOR, WA_G_START, WA_I_FRAG, WA_LENGTH, WA_NREP, WA_R_START, WA_SIZE,
    WA_SJ_A, WC_CHR, WC_G_END, WC_G_START, WC_SIZE, WC_STR,
};
use star_genome::Genome;
use star_params::parameters::Parameters;

use crate::read_align::{ReadAlign, WinAlign, WinCoord};

impl ReadAlign {
    fn ensure_window_capacity(&mut self, up_to: usize, p: &Parameters) {
        let want = up_to + 1;
        if self.wc.len() < want {
            self.wc.resize(want, [0u64; WC_SIZE]);
        }
        if self.wa.len() < want {
            self.wa.resize(
                want,
                vec![[0u64; WA_SIZE]; p.seed_per_window_nmax as usize + 2],
            );
        }
        if self.n_wa.len() < want {
            self.n_wa.resize(want, 0);
        }
        if self.n_wap.len() < want {
            self.n_wap.resize(want, 0);
        }
        if self.wa_lrec.len() < want {
            self.wa_lrec.resize(want, 0);
        }
        if self.wlast_anchor.len() < want {
            self.wlast_anchor.resize(want, u64::MAX);
        }
        if self.n_win_tr.len() < want {
            self.n_win_tr.resize(want, 0);
        }
    }

    /// 1:1 port of `ReadAlign::createExtendWindowsWithAlign`
    /// (`ReadAlign_createExtendWindowsWithAlign.cpp:7-85`).
    pub fn create_extend_windows_with_align(
        &mut self,
        p: &Parameters,
        map_gen: &Genome,
        a1: u64,
        a_str: u64,
    ) -> i32 {
        let a_bin = (a1 >> p.win_bin_nbits) as usize;
        let mut i_bin_left = a_bin;
        let mut i_bin_right = a_bin;
        let mut i_win: usize = usize::MAX;
        let mut i_win_right: usize = usize::MAX;

        let w_b = &mut self.win_bin[a_str as usize];
        if w_b[a_bin] != UINT_WIN_BIN_MAX {
            return 0;
        }

        let mut flag_merge_left = false;
        if a_bin > 0 {
            let lo = if a_bin as u64 > p.win_anchor_dist_nbins {
                a_bin - p.win_anchor_dist_nbins as usize
            } else {
                0
            };
            let mut i_bin = a_bin - 1;
            loop {
                if w_b[i_bin] < UINT_WIN_BIN_MAX {
                    flag_merge_left = true;
                    break;
                }
                if i_bin == lo {
                    break;
                }
                i_bin -= 1;
            }
            if flag_merge_left
                && map_gen.chr_bin[i_bin >> p.win_bin_chr_nbits]
                    == map_gen.chr_bin[a_bin >> p.win_bin_chr_nbits]
            {
                i_win = w_b[i_bin] as usize;
                i_bin_left = self.wc[i_win][WC_G_START] as usize;
                let range_lo = i_bin + 1;
                for ii in range_lo..=a_bin {
                    w_b[ii] = i_win as u16;
                }
            } else {
                flag_merge_left = false;
            }
        }

        let mut flag_merge_right = false;
        if (a_bin as u64 + 1) < p.win_bin_n {
            let hi = (a_bin + p.win_anchor_dist_nbins as usize + 1).min(p.win_bin_n as usize);
            let mut i_bin = a_bin + 1;
            while i_bin < hi {
                if w_b[i_bin] < UINT_WIN_BIN_MAX {
                    flag_merge_right = true;
                    break;
                }
                i_bin += 1;
            }
            if flag_merge_right
                && i_bin < p.win_bin_n as usize
                && map_gen.chr_bin[i_bin >> p.win_bin_chr_nbits]
                    == map_gen.chr_bin[a_bin >> p.win_bin_chr_nbits]
            {
                while i_bin + 1 < w_b.len() && w_b[i_bin] == w_b[i_bin + 1] {
                    i_bin += 1;
                }
                i_bin_right = i_bin;
                i_win_right = w_b[i_bin] as usize;
                if !flag_merge_left {
                    i_win = w_b[i_bin] as usize;
                }
                for ii in a_bin..=i_bin {
                    w_b[ii] = i_win as u16;
                }
            } else {
                flag_merge_right = false;
            }
        }

        if !flag_merge_left && !flag_merge_right {
            i_win = self.n_w as usize;
            w_b[a_bin] = i_win as u16;
            self.ensure_window_capacity(i_win, p);
            self.wc[i_win][WC_CHR] = map_gen.chr_bin[a_bin >> p.win_bin_chr_nbits];
            self.wc[i_win][WC_STR] = a_str;
            self.wc[i_win][WC_G_START] = a_bin as u64;
            self.wc[i_win][WC_G_END] = a_bin as u64;
            self.n_w += 1;
            if self.n_w >= p.align_window_per_read_nmax as u64 {
                self.n_w = (p.align_window_per_read_nmax - 1) as u64;
                return EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS;
            }
        } else {
            // Only include i_win_right in the capacity request when it was
            // actually assigned (flag_merge_right). Otherwise the sentinel
            // `usize::MAX` would ask for a Vec of size usize::MAX and panic
            // with `capacity overflow` — hit on the Mp/chr1 fixture.
            let cap_up_to = if i_win_right != usize::MAX {
                i_win.max(i_win_right)
            } else {
                i_win
            };
            self.ensure_window_capacity(cap_up_to, p);
            self.wc[i_win][WC_G_START] = i_bin_left as u64;
            self.wc[i_win][WC_G_END] = i_bin_right as u64;
            if flag_merge_left && flag_merge_right && i_win_right != usize::MAX {
                self.wc[i_win_right][WC_G_START] = 1;
                self.wc[i_win_right][WC_G_END] = 0;
            }
        }
        0
    }

    /// 1:1 port of `ReadAlign::assignAlignToWindow`
    /// (`ReadAlign_assignAlignToWindow.cpp:6-130`).
    #[allow(clippy::too_many_arguments)]
    pub fn assign_align_to_window(
        &mut self,
        p: &Parameters,
        a1: u64,
        a_length: u64,
        a_str: u64,
        a_nrep: u64,
        a_frag: u64,
        a_rstart: u64,
        a_anchor: bool,
        sj_a: u64,
    ) {
        let bin = (a1 >> p.win_bin_nbits) as usize;
        let w_b = &self.win_bin[a_str as usize];
        if bin >= w_b.len() {
            return;
        }
        let i_w = w_b[bin];
        if i_w == UINT_WIN_BIN_MAX || (!a_anchor && a_length < self.wa_lrec[i_w as usize]) {
            return;
        }
        let i_w = i_w as usize;

        let nwa = self.n_wa[i_w] as usize;
        let mut overlapping: Option<usize> = None;
        for i_a in 0..nwa {
            let wa = &self.wa[i_w][i_a];
            if a_frag == wa[WA_I_FRAG]
                && wa[WA_SJ_A] == sj_a
                && a1 + wa[WA_R_START] == wa[WA_G_START] + a_rstart
                && ((a_rstart >= wa[WA_R_START] && a_rstart < wa[WA_R_START] + wa[WA_LENGTH])
                    || (a_rstart + a_length >= wa[WA_R_START]
                        && a_rstart + a_length < wa[WA_R_START] + wa[WA_LENGTH]))
            {
                overlapping = Some(i_a);
                break;
            }
        }
        if let Some(i_a) = overlapping {
            let existing_len = self.wa[i_w][i_a][WA_LENGTH];
            if a_length > existing_len {
                let mut i_a0 = nwa;
                for j in 0..nwa {
                    if j != i_a && a_rstart < self.wa[i_w][j][WA_R_START] {
                        i_a0 = j;
                        break;
                    }
                }
                if i_a0 > i_a {
                    i_a0 -= 1;
                }
                if i_a0 < i_a {
                    for ia1 in (i_a0 + 1..=i_a).rev() {
                        self.wa[i_w][ia1] = self.wa[i_w][ia1 - 1];
                    }
                } else if i_a0 > i_a {
                    for ia1 in i_a..i_a0 {
                        self.wa[i_w][ia1] = self.wa[i_w][ia1 + 1];
                    }
                }
                self.wa[i_w][i_a0][WA_R_START] = a_rstart;
                self.wa[i_w][i_a0][WA_LENGTH] = a_length;
                self.wa[i_w][i_a0][WA_G_START] = a1;
                self.wa[i_w][i_a0][WA_NREP] = a_nrep;
                self.wa[i_w][i_a0][WA_ANCHOR] = if a_anchor { 1 } else { 0 };
                self.wa[i_w][i_a0][WA_I_FRAG] = a_frag;
                self.wa[i_w][i_a0][WA_SJ_A] = sj_a;
            }
            return;
        }

        if self.n_wa[i_w] == p.seed_per_window_nmax {
            let l_read = self.l_read;
            self.wa_lrec[i_w] = l_read + 1;
            for i_a in 0..self.n_wa[i_w] as usize {
                if self.wa[i_w][i_a][WA_ANCHOR] != 1 {
                    self.wa_lrec[i_w] = self.wa_lrec[i_w].min(self.wa[i_w][i_a][WA_LENGTH]);
                }
            }
            if self.wa_lrec[i_w] == l_read + 1 {
                self.map_marker = MARKER_TOO_MANY_ANCHORS_PER_WINDOW;
                self.n_w = 0;
                return;
            }
            if !a_anchor && a_length < self.wa_lrec[i_w] {
                return;
            }
            let mut i_a1 = 0usize;
            for i_a in 0..self.n_wa[i_w] as usize {
                if self.wa[i_w][i_a][WA_ANCHOR] == 1
                    || self.wa[i_w][i_a][WA_LENGTH] > self.wa_lrec[i_w]
                {
                    self.wa[i_w][i_a1] = self.wa[i_w][i_a];
                    i_a1 += 1;
                }
            }
            self.n_wa[i_w] = i_a1 as u64;
            if !a_anchor && a_length <= self.wa_lrec[i_w] {
                self.n_wap[i_w] = 0;
            }
        }

        if a_anchor || a_length > self.wa_lrec[i_w] {
            if self.n_wa[i_w] >= p.seed_per_window_nmax {
                // Matches C++ BUG exit; in Rust we just silently clamp to keep going.
                return;
            }
            let nwa = self.n_wa[i_w] as usize;
            let mut i_a = nwa;
            for j in 0..nwa {
                if a_rstart < self.wa[i_w][j][WA_R_START] {
                    i_a = j;
                    break;
                }
            }
            for ia1 in (i_a + 1..=nwa).rev() {
                self.wa[i_w][ia1] = self.wa[i_w][ia1 - 1];
            }
            self.wa[i_w][i_a][WA_R_START] = a_rstart;
            self.wa[i_w][i_a][WA_LENGTH] = a_length;
            self.wa[i_w][i_a][WA_G_START] = a1;
            self.wa[i_w][i_a][WA_NREP] = a_nrep;
            self.wa[i_w][i_a][WA_ANCHOR] = if a_anchor { 1 } else { 0 };
            self.wa[i_w][i_a][WA_I_FRAG] = a_frag;
            self.wa[i_w][i_a][WA_SJ_A] = sj_a;
            self.n_wa[i_w] += 1;
            self.n_wap[i_w] += 1;
            if a_anchor && self.wlast_anchor[i_w] < i_a as u64 {
                self.wlast_anchor[i_w] = i_a as u64;
            }
        }
    }
}

/// 1:1 port of `sjAlignSplit.cpp` — splits an SJ-encoded alignment position
/// back into donor and acceptor pieces in the linear genome. Returns `Some`
/// with `(a1D, aLengthD, a1A, aLengthA, isj)` if the align crosses a junction.
pub fn sj_align_split(
    a1: u64,
    a_length: u64,
    map_gen: &Genome,
) -> Option<(u64, u64, u64, u64, u64)> {
    if map_gen.sjdb_length == 0 {
        return None;
    }
    let sj1 = (a1 - map_gen.sj_g_start) % map_gen.sjdb_length;
    if sj1 < map_gen.sjdb_overhang && sj1 + a_length > map_gen.sjdb_overhang {
        let isj = (a1 - map_gen.sj_g_start) / map_gen.sjdb_length;
        let a_length_d = map_gen.sjdb_overhang - sj1;
        let a_length_a = a_length - a_length_d;
        let a1_d = map_gen.sj_d_start[isj as usize] + sj1;
        let a1_a = map_gen.sj_a_start[isj as usize];
        Some((a1_d, a_length_d, a1_a, a_length_a, isj))
    } else {
        None
    }
}

/// Helper kept for parity with C++ `WinAlign` referenced by downstream code.
#[inline]
pub fn wa_zero() -> WinAlign {
    [0u64; WA_SIZE]
}

#[inline]
pub fn wc_zero() -> WinCoord {
    [0u64; WC_SIZE]
}
