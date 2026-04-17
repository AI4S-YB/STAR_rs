//! Stitch pipeline — ports the following C++ sources into a single module:
//!
//! - `extendAlign.cpp`             -> [`extend_align`]
//! - `binarySearch2.cpp`           -> [`binary_search2`]
//! - `blocksOverlap.cpp`           -> [`blocks_overlap`]
//! - `stitchAlignToTranscript.cpp` -> [`stitch_align_to_transcript`]
//!
//! All functions are pure; they operate on raw genome / read pointers plus
//! the [`Transcript`] state. Unsafe blocks are kept narrow and mirror the C++
//! pointer arithmetic 1:1.

use star_core::types::{
    DEF_READ_SEQ_LENGTH_MAX, EX_G, EX_I_FRAG, EX_L, EX_R, EX_SJ_A, MARK_FRAG_SPACER_BASE,
    MAX_N_EXONS, MAX_SJ_REPEAT_SEARCH, SCORE_MATCH,
};
use star_core::UInt;
use star_genome::Genome;
use star_params::parameters::Parameters;

use crate::transcript::Transcript;

/// 1:1 port of `binarySearch2.cpp`. Returns the index of the matching (x, y)
/// pair, `-1` if not present, or `-2` for overflow (should not happen in
/// well-formed data).
///
/// Note: parameters use `u64` for (x, y) and `i64` for the returned index to
/// stay in sync with C++ `uint` (`unsigned long long`) / `int`.
pub fn binary_search2(x: u64, y: u64, xs: &[u64], ys: &[u64]) -> i64 {
    let n = xs.len();
    if n == 0 || x > xs[n - 1] || x < xs[0] {
        return -1;
    }
    let mut i1: usize = 0;
    let mut i2: usize = n - 1;
    let mut i3: usize;
    while i2 > i1 + 1 {
        i3 = (i1 + i2) / 2;
        if xs[i3] > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    if x == xs[i1] {
        i3 = i1;
    } else if x == xs[i2] {
        i3 = i2;
    } else {
        return -1;
    }
    let mut jj: isize = i3 as isize;
    while jj >= 0 {
        if x != xs[jj as usize] {
            break;
        }
        if y == ys[jj as usize] {
            return jj as i64;
        }
        jj -= 1;
    }
    let mut jj = i3;
    while jj < n {
        if x != xs[jj] {
            return -1;
        }
        if y == ys[jj] {
            return jj as i64;
        }
        jj += 1;
    }
    -2
}

/// 1:1 port of `blocksOverlap.cpp` — overlap of two transcripts' exon blocks.
pub fn blocks_overlap(t1: &Transcript, t2: &Transcript) -> u64 {
    let mut i1: usize = 0;
    let mut i2: usize = 0;
    let mut n_overlap: u64 = 0;
    while (i1 as u64) < t1.n_exons && (i2 as u64) < t2.n_exons {
        let rs1 = t1.exons[i1][EX_R];
        let rs2 = t2.exons[i2][EX_R];
        let re1 = t1.exons[i1][EX_R] + t1.exons[i1][EX_L];
        let re2 = t2.exons[i2][EX_R] + t2.exons[i2][EX_L];
        let gs1 = t1.exons[i1][EX_G];
        let gs2 = t2.exons[i2][EX_G];

        if rs1 >= re2 {
            i2 += 1;
        } else if rs2 >= re1 {
            i1 += 1;
        } else if gs1.wrapping_sub(rs1) != gs2.wrapping_sub(rs2) {
            if re1 >= re2 {
                i2 += 1;
            }
            if re2 >= re1 {
                i1 += 1;
            }
        } else {
            n_overlap += re1.min(re2) - rs1.max(rs2);
            if re1 >= re2 {
                i2 += 1;
            }
            if re2 >= re1 {
                i1 += 1;
            }
        }
    }
    n_overlap
}

/// 1:1 port of `extendAlign.cpp`.
///
/// - `r` / `g` are raw base pointers; negative offsets are supported through
///   the pre-padded genome buffer (see `LOAD_L` in `star_genome::load`).
/// - Returns `true` if any extension landed (matches C++ `extDone`).
///
/// # Safety
/// - `r` must be valid for reads of `[r_start, r_start + L * |dR|)`.
/// - `g` must be valid for reads of `[g_start + min(0, dG*L), g_start + max(0, dG*L))`
///   within the loaded genome buffer including the `LOAD_L` padding.
#[allow(clippy::too_many_arguments)]
pub unsafe fn extend_align(
    r: *const u8,
    g: *const u8,
    r_start: u64,
    g_start: u64,
    d_r: i64,
    d_g: i64,
    l: u64,
    l_prev: u64,
    n_mm_prev: u64,
    n_mm_max: u64,
    p_mm_max: f64,
    extend_to_end: bool,
    tr_a: &mut Transcript,
) -> bool {
    let mut score: i64 = 0;
    let mut n_match: u64 = 0;
    let mut n_mm: u64 = 0;
    tr_a.max_score = 0;

    let r = r.add(r_start as usize);
    let g = g.add(g_start as usize);

    if extend_to_end {
        let mut i_ext: i64 = 0;
        while i_ext < l as i64 {
            let i_s = d_r * i_ext;
            let i_g = d_g * i_ext;
            let g_byte = *g.offset(i_g as isize);
            if (g_start as i64 + i_g) == -1_i64 || g_byte == 5 {
                tr_a.extend_l = 0;
                tr_a.max_score = -999_999_999;
                tr_a.n_match = 0;
                tr_a.n_mm = n_mm_max + 1;
                return true;
            }
            let r_byte = *r.offset(i_s as isize);
            if r_byte == MARK_FRAG_SPACER_BASE {
                break;
            }
            if r_byte > 3 || g_byte > 3 {
                i_ext += 1;
                continue;
            }
            if g_byte == r_byte {
                n_match += 1;
                score += SCORE_MATCH as i64;
            } else {
                n_mm += 1;
                score -= SCORE_MATCH as i64;
            }
            i_ext += 1;
        }
        return if i_ext > 0 {
            tr_a.extend_l = i_ext as u64;
            tr_a.max_score = score as i32;
            tr_a.n_match = n_match;
            tr_a.n_mm = n_mm;
            true
        } else {
            false
        };
    }

    for i in 0..l as i64 {
        let i_s = d_r * i;
        let i_g = d_g * i;
        let g_byte = *g.offset(i_g as isize);
        let r_byte = *r.offset(i_s as isize);
        if (g_start as i64 + i_g) == -1_i64 || g_byte == 5 || r_byte == MARK_FRAG_SPACER_BASE {
            break;
        }
        if r_byte > 3 || g_byte > 3 {
            continue;
        }
        if g_byte == r_byte {
            n_match += 1;
            score += SCORE_MATCH as i64;
            if score > tr_a.max_score as i64 {
                let limit =
                    (p_mm_max * (l_prev as i64 + i + 1) as f64).min(n_mm_max as f64) as u64;
                if n_mm + n_mm_prev <= limit {
                    tr_a.extend_l = (i + 1) as u64;
                    tr_a.max_score = score as i32;
                    tr_a.n_match = n_match;
                    tr_a.n_mm = n_mm;
                }
            }
        } else {
            let limit = (p_mm_max * (l_prev + l) as f64).min(n_mm_max as f64) as u64;
            if n_mm + n_mm_prev >= limit {
                break;
            }
            n_mm += 1;
            score -= SCORE_MATCH as i64;
        }
    }

    tr_a.extend_l > 0
}

/// 1:1 port of `stitchAlignToTranscript.cpp`.
///
/// # Safety
/// `r` and `g` must point into properly allocated read/genome buffers
/// such that all offsets computed below (including those which wrap into the
/// pre-load `LOAD_L` pad) are valid.
#[allow(clippy::too_many_arguments)]
pub unsafe fn stitch_align_to_transcript(
    r_aend: u64,
    g_aend: u64,
    mut r_bstart: u64,
    mut g_bstart: u64,
    mut l: u64,
    i_frag_b: u64,
    sj_ab: u64,
    p: &Parameters,
    r: *const u8,
    map_gen: &Genome,
    tr_a: &mut Transcript,
    out_filter_mismatch_nmax_total: u64,
) -> i32 {
    if tr_a.n_exons as usize >= MAX_N_EXONS {
        return -1_000_010;
    }
    let g = map_gen.g.as_ptr().add(star_genome::load::LOAD_L);
    let mut score: i64 = 0;
    let scored: i64 = SCORE_MATCH as i64;

    let sj_none: u64 = u64::MAX;

    let nexons = tr_a.n_exons as usize;
    if sj_ab != sj_none
        && tr_a.exons[nexons - 1][EX_SJ_A] == sj_ab
        && tr_a.exons[nexons - 1][EX_I_FRAG] == i_frag_b
        && r_bstart == r_aend + 1
        && g_aend + 1 < g_bstart
    {
        let sj_idx = sj_ab as usize;
        if map_gen.sjdb_motif[sj_idx] == 0
            && (l <= map_gen.sjdb_shift_right[sj_idx] as u64
                || tr_a.exons[nexons - 1][EX_L] <= map_gen.sjdb_shift_left[sj_idx] as u64)
        {
            return -1_000_006;
        }
        tr_a.exons[nexons][EX_L] = l;
        tr_a.exons[nexons][EX_R] = r_bstart;
        tr_a.exons[nexons][EX_G] = g_bstart;
        tr_a.canon_sj[nexons - 1] = map_gen.sjdb_motif[sj_idx] as i32;
        tr_a.shift_sj[nexons - 1][0] = map_gen.sjdb_shift_left[sj_idx] as UInt;
        tr_a.shift_sj[nexons - 1][1] = map_gen.sjdb_shift_right[sj_idx] as UInt;
        tr_a.sj_annot[nexons - 1] = 1;
        tr_a.sj_str[nexons - 1] = map_gen.sjdb_strand[sj_idx];
        tr_a.n_exons += 1;
        tr_a.n_match += l;
        score += (l as i64) * scored;
        score += p.p_ge.sjdb_score as i64;
    } else {
        tr_a.sj_annot[nexons - 1] = 0;
        tr_a.sj_str[nexons - 1] = 0;

        if tr_a.exons[nexons - 1][EX_I_FRAG] == i_frag_b {
            let g_bend = g_bstart + l - 1;
            let r_bend = r_bstart + l - 1;

            if r_bend <= r_aend {
                return -1_000_001;
            }
            if g_bend <= g_aend && tr_a.exons[nexons - 1][EX_I_FRAG] == i_frag_b {
                return -1_000_002;
            }

            if r_bstart <= r_aend {
                g_bstart += r_aend - r_bstart + 1;
                r_bstart = r_aend + 1;
                l = r_bend - r_bstart + 1;
            }

            score += ((r_bend - r_bstart + 1) as i64) * scored;

            let g_gap: i64 = g_bstart as i64 - g_aend as i64 - 1;
            let r_gap: i64 = r_bstart as i64 - r_aend as i64 - 1;

            let mut n_match: u64 = l;
            let mut n_mm: u64 = 0;
            let mut del_l: i64 = 0;
            let mut ins_l: i64 = 0;
            let mut n_ins: u64 = 0;
            let mut n_del: u64 = 0;
            let mut j_r: i64 = 0;
            let mut j_can: i32 = 999;
            let g_bstart1: i64 = g_bstart as i64 - r_gap - 1;

            if g_gap == 0 && r_gap == 0 {
                // no-op
            } else if g_gap > 0 && r_gap > 0 && r_gap == g_gap {
                for ii in 1..=r_gap {
                    let g_byte = *g.offset(g_aend as isize + ii as isize);
                    let r_byte = *r.offset(r_aend as isize + ii as isize);
                    if g_byte < 4 && r_byte < 4 {
                        if r_byte == g_byte {
                            score += scored;
                            n_match += 1;
                        } else {
                            score -= scored;
                            n_mm += 1;
                        }
                    }
                }
            } else if g_gap > r_gap {
                n_del = 1;
                del_l = g_gap - r_gap;
                if del_l as u64 > p.align_intron_max && p.align_intron_max > 0 {
                    return -1_000_003;
                }

                let mut score1: i64 = 0;
                let mut j_r1: i64 = 1;
                loop {
                    j_r1 -= 1;
                    let r_byte = *r.offset(r_aend as isize + j_r1 as isize);
                    let g1_byte = *g.offset(g_bstart1 as isize + j_r1 as isize);
                    let g2_byte = *g.offset(g_aend as isize + j_r1 as isize);
                    if r_byte != g1_byte && g1_byte < 4 && r_byte == g2_byte {
                        score1 -= scored;
                    }
                    let keep = score1 + p.score_stitch_sj_shift as i64 >= 0
                        && tr_a.exons[nexons - 1][EX_L] as i64 + j_r1 > 1;
                    if !keep {
                        break;
                    }
                }

                let mut max_score2: i64 = -999_999;
                let mut score1: i64 = 0;
                let mut j_pen: i32 = 0;
                loop {
                    let r_byte = *r.offset(r_aend as isize + j_r1 as isize);
                    let g_a = *g.offset(g_aend as isize + j_r1 as isize);
                    let g_b = *g.offset(g_bstart1 as isize + j_r1 as isize);
                    if r_byte == g_a && r_byte != g_b {
                        score1 += scored;
                    }
                    if r_byte != g_a && r_byte == g_b {
                        score1 -= scored;
                    }

                    let mut j_can1: i32 = -1;
                    let mut j_pen1: i32 = 0;
                    let mut score2 = score1;

                    if del_l as u64 >= p.align_intron_min {
                        let g_a1 = *g.offset(g_aend as isize + j_r1 as isize + 1);
                        let g_a2 = *g.offset(g_aend as isize + j_r1 as isize + 2);
                        let g_bm1 = *g.offset(g_bstart1 as isize + j_r1 as isize - 1);
                        let g_b0 = *g.offset(g_bstart1 as isize + j_r1 as isize);

                        if g_a1 == 2 && g_a2 == 3 && g_bm1 == 0 && g_b0 == 2 {
                            j_can1 = 1;
                        } else if g_a1 == 1 && g_a2 == 3 && g_bm1 == 0 && g_b0 == 1 {
                            j_can1 = 2;
                        } else if g_a1 == 2 && g_a2 == 1 && g_bm1 == 0 && g_b0 == 2 {
                            j_can1 = 3;
                            j_pen1 = p.score_gap_gcag;
                        } else if g_a1 == 1 && g_a2 == 3 && g_bm1 == 2 && g_b0 == 1 {
                            j_can1 = 4;
                            j_pen1 = p.score_gap_gcag;
                        } else if g_a1 == 0 && g_a2 == 3 && g_bm1 == 0 && g_b0 == 1 {
                            j_can1 = 5;
                            j_pen1 = p.score_gap_atac;
                        } else if g_a1 == 2 && g_a2 == 3 && g_bm1 == 0 && g_b0 == 3 {
                            j_can1 = 6;
                            j_pen1 = p.score_gap_atac;
                        } else {
                            j_can1 = 0;
                            j_pen1 = p.score_gap_noncan;
                        }
                        score2 += j_pen1 as i64;
                    }

                    if max_score2 < score2 {
                        max_score2 = score2;
                        j_r = j_r1;
                        j_can = j_can1;
                        j_pen = j_pen1;
                    }
                    j_r1 += 1;
                    if j_r1 >= r_bend as i64 - r_aend as i64 {
                        break;
                    }
                }

                let mut jj_l: u64 = 0;
                let mut jj_r: u64 = 0;
                while g_aend as i64 + j_r >= jj_l as i64
                    && *g.offset(g_aend as isize - jj_l as isize + j_r as isize)
                        == *g.offset(g_bstart1 as isize - jj_l as isize + j_r as isize)
                    && *g.offset(g_aend as isize - jj_l as isize + j_r as isize) < 4
                    && (jj_l as usize) <= MAX_SJ_REPEAT_SEARCH
                {
                    jj_l += 1;
                }
                while (g_aend as i64 + jj_r as i64 + j_r + 1) < map_gen.n_genome as i64
                    && *g.offset(g_aend as isize + jj_r as isize + j_r as isize + 1)
                        == *g.offset(g_bstart1 as isize + jj_r as isize + j_r as isize + 1)
                    && *g.offset(g_aend as isize + jj_r as isize + j_r as isize + 1) < 4
                    && (jj_r as usize) <= MAX_SJ_REPEAT_SEARCH
                {
                    jj_r += 1;
                }

                if j_can <= 0 {
                    j_r -= jj_l as i64;
                    if tr_a.exons[nexons - 1][EX_L] as i64 + j_r < 1 {
                        return -1_000_005;
                    }
                    jj_r += jj_l;
                    jj_l = 0;
                }

                let ii_start = 1_i64.min(j_r + 1);
                let ii_end = r_gap.max(j_r);
                for ii in ii_start..=ii_end {
                    let g1: i64 = if ii <= j_r {
                        g_aend as i64 + ii
                    } else {
                        g_bstart1 + ii
                    };
                    let g_byte = *g.offset(g1 as isize);
                    let r_byte = *r.offset(r_aend as isize + ii as isize);
                    if g_byte < 4 && r_byte < 4 {
                        if r_byte == g_byte {
                            if ii >= 1 && ii <= r_gap {
                                score += scored;
                                n_match += 1;
                            }
                        } else {
                            score -= scored;
                            n_mm += 1;
                            if ii < 1 || ii > r_gap {
                                score -= scored;
                                n_match = n_match.saturating_sub(1);
                            }
                        }
                    }
                }

                if map_gen.sjdb_n > 0 {
                    let j_s = (g_aend as i64 + j_r + 1) as u64;
                    let j_e = (g_bstart1 + j_r) as u64;
                    let sjdb_ind =
                        binary_search2(j_s, j_e, &map_gen.sjdb_start, &map_gen.sjdb_end);
                    if sjdb_ind < 0 {
                        if del_l as u64 >= p.align_intron_min {
                            score += (p.score_gap + j_pen) as i64;
                        } else {
                            score += del_l * p.score_del_base as i64 + p.score_del_open as i64;
                            j_can = -1;
                            tr_a.sj_annot[nexons - 1] = 0;
                        }
                    } else {
                        let si = sjdb_ind as usize;
                        j_can = map_gen.sjdb_motif[si] as i32;
                        if map_gen.sjdb_motif[si] == 0 {
                            if l <= map_gen.sjdb_shift_left[si] as u64
                                || tr_a.exons[nexons - 1][EX_L]
                                    <= map_gen.sjdb_shift_left[si] as u64
                            {
                                return -1_000_006;
                            }
                            j_r += map_gen.sjdb_shift_left[si] as i64;
                            if r_aend as i64 + j_r >= r_bend as i64 {
                                return -1_000_006;
                            }
                            jj_l = map_gen.sjdb_shift_left[si] as u64;
                            jj_r = map_gen.sjdb_shift_right[si] as u64;
                        }
                        tr_a.sj_annot[nexons - 1] = 1;
                        tr_a.sj_str[nexons - 1] = map_gen.sjdb_strand[si];
                        score += p.p_ge.sjdb_score as i64;
                    }
                } else if del_l as u64 >= p.align_intron_min {
                    score += (p.score_gap + j_pen) as i64;
                } else {
                    score += del_l * p.score_del_base as i64 + p.score_del_open as i64;
                    j_can = -1;
                    tr_a.sj_annot[nexons - 1] = 0;
                }

                tr_a.shift_sj[nexons - 1][0] = jj_l;
                tr_a.shift_sj[nexons - 1][1] = jj_r;
                tr_a.canon_sj[nexons - 1] = j_can;

                if tr_a.sj_annot[nexons - 1] == 0 {
                    tr_a.sj_str[nexons - 1] = if j_can > 0 {
                        (2 - (j_can % 2)) as u8
                    } else {
                        0
                    };
                }
            } else if r_gap > g_gap {
                ins_l = r_gap - g_gap;
                n_ins = 1;
                if g_gap == 0 {
                    j_r = 0;
                } else if g_gap < 0 {
                    j_r = 0;
                    for _ii in 0..(-g_gap) {
                        score -= scored;
                    }
                } else {
                    let mut score1: i64 = 0;
                    let mut max_score1: i64 = 0;
                    for j_r1 in 1..=g_gap {
                        let g_byte = *g.offset(g_aend as isize + j_r1 as isize);
                        if g_byte < 4 {
                            let r0 = *r.offset(r_aend as isize + j_r1 as isize);
                            let r1 = *r.offset(r_aend as isize + ins_l as isize + j_r1 as isize);
                            if r0 == g_byte {
                                score1 += scored;
                            } else {
                                score1 -= scored;
                            }
                            if r1 == g_byte {
                                score1 -= scored;
                            } else {
                                score1 += scored;
                            }
                        }
                        if score1 > max_score1
                            || (score1 == max_score1 && p.align_insertion_flush.flush_right)
                        {
                            max_score1 = score1;
                            j_r = j_r1;
                        }
                    }
                    for ii in 1..=g_gap {
                        let r1 = r_aend_plus_ins(r_aend, ii, j_r, ins_l);
                        let g_byte = *g.offset(g_aend as isize + ii as isize);
                        let r_byte = *r.offset(r1 as isize);
                        if g_byte < 4 && r_byte < 4 {
                            if r_byte == g_byte {
                                score += scored;
                                n_match += 1;
                            } else {
                                score -= scored;
                                n_mm += 1;
                            }
                        }
                    }
                }

                if p.align_insertion_flush.flush_right {
                    while j_r < r_bend as i64 - r_aend as i64 - ins_l {
                        let r_byte = *r.offset(r_aend as isize + j_r as isize + 1);
                        let g_byte = *g.offset(g_aend as isize + j_r as isize + 1);
                        if r_byte != g_byte || g_byte == 4 {
                            break;
                        }
                        j_r += 1;
                    }
                    if j_r == r_bend as i64 - r_aend as i64 - ins_l {
                        return -1_000_009;
                    }
                }
                score += ins_l * p.score_ins_base as i64 + p.score_ins_open as i64;
                j_can = -2;
            }

            // Mirror C++ `#ifdef COMPILE_FOR_LONG_READS` — short reads path
            // additionally rejects excessive MM for non-GT/AG junctions.
            let ok = {
                let idx = ((j_can + 1) / 2).max(0) as usize;
                // C++ casts `P.alignSJstitchMismatchNmax[idx]` to `uint`, so
                // the sentinel `-1` becomes `UINT_MAX` and the inequality is
                // effectively "no limit". Preserve that semantics here.
                let ok_mm = if idx < p.align_sj_stitch_mismatch_nmax.len() {
                    let limit = p.align_sj_stitch_mismatch_nmax[idx];
                    limit < 0 || n_mm as i64 <= limit
                } else {
                    false
                };
                tr_a.n_mm + n_mm <= out_filter_mismatch_nmax_total
                    && (j_can < 0 || (j_can < 7 && ok_mm))
            };

            if ok {
                tr_a.n_mm += n_mm;
                tr_a.n_match += n_match;
                if del_l as u64 >= p.align_intron_min {
                    tr_a.n_gap += n_del;
                    tr_a.l_gap += del_l as u64;
                } else {
                    tr_a.n_del += n_del;
                    tr_a.l_del += del_l as u64;
                }

                if del_l == 0 && ins_l == 0 {
                    tr_a.exons[nexons - 1][EX_L] += r_bend - r_aend;
                } else if del_l > 0 {
                    tr_a.exons[nexons - 1][EX_L] = (tr_a.exons[nexons - 1][EX_L] as i64 + j_r) as u64;
                    tr_a.exons[nexons][EX_L] = (r_bend as i64 - r_aend as i64 - j_r) as u64;
                    tr_a.exons[nexons][EX_R] = (r_aend as i64 + j_r + 1) as u64;
                    tr_a.exons[nexons][EX_G] = (g_bstart1 + j_r + 1) as u64;
                    tr_a.n_exons += 1;
                } else if ins_l > 0 {
                    tr_a.n_ins += n_ins;
                    tr_a.l_ins += ins_l as u64;
                    tr_a.exons[nexons - 1][EX_L] = (tr_a.exons[nexons - 1][EX_L] as i64 + j_r) as u64;
                    tr_a.exons[nexons][EX_L] = (r_bend as i64 - r_aend as i64 - j_r - ins_l) as u64;
                    tr_a.exons[nexons][EX_R] = (r_aend as i64 + j_r + ins_l + 1) as u64;
                    tr_a.exons[nexons][EX_G] = g_aend + 1 + j_r as u64;
                    tr_a.canon_sj[nexons - 1] = -2;
                    tr_a.sj_annot[nexons - 1] = 0;
                    tr_a.n_exons += 1;
                }
            } else {
                return -1_000_007;
            }

            let _ = g_bend; // kept for parity with C++ symbol usage
        } else if g_bstart + tr_a.exons[0][EX_R] + p.align_ends_protrude.n_bases_max as u64
            >= tr_a.exons[0][EX_G]
            || tr_a.exons[0][EX_G] < tr_a.exons[0][EX_R]
        {
            if p.align_mates_gap_max > 0
                && g_bstart
                    > tr_a.exons[nexons - 1][EX_G]
                        + tr_a.exons[nexons - 1][EX_L]
                        + p.align_mates_gap_max
            {
                return -1_000_004;
            }

            score += (l as i64) * scored;

            let mut tr_extend = Transcript::new();
            tr_extend.reset();
            if extend_align(
                r,
                g,
                r_aend + 1,
                g_aend + 1,
                1,
                1,
                DEF_READ_SEQ_LENGTH_MAX as u64,
                tr_a.n_match,
                tr_a.n_mm,
                out_filter_mismatch_nmax_total,
                p.out_filter_mismatch_nover_lmax,
                p.align_ends_type.ext[tr_a.exons[nexons - 1][EX_I_FRAG] as usize][1],
                &mut tr_extend,
            ) {
                tr_a.add(&tr_extend);
                score += tr_extend.max_score as i64;
                tr_a.exons[nexons - 1][EX_L] += tr_extend.extend_l;
            }

            tr_a.exons[nexons][EX_R] = r_bstart;
            tr_a.exons[nexons][EX_G] = g_bstart;
            tr_a.exons[nexons][EX_L] = l;
            tr_a.n_match += l;

            tr_extend.reset();
            let extlen: u64 = if p.align_ends_type.ext[i_frag_b as usize][1] {
                DEF_READ_SEQ_LENGTH_MAX as u64
            } else {
                g_bstart - tr_a.exons[0][EX_G] + tr_a.exons[0][EX_R]
            };
            if extend_align(
                r,
                g,
                r_bstart - 1,
                g_bstart - 1,
                -1,
                -1,
                extlen,
                tr_a.n_match,
                tr_a.n_mm,
                out_filter_mismatch_nmax_total,
                p.out_filter_mismatch_nover_lmax,
                p.align_ends_type.ext[i_frag_b as usize][1],
                &mut tr_extend,
            ) {
                tr_a.add(&tr_extend);
                score += tr_extend.max_score as i64;
                tr_a.exons[nexons][EX_R] -= tr_extend.extend_l;
                tr_a.exons[nexons][EX_G] -= tr_extend.extend_l;
                tr_a.exons[nexons][EX_L] += tr_extend.extend_l;
            }

            tr_a.canon_sj[nexons - 1] = -3;
            tr_a.sj_annot[nexons - 1] = 0;
            tr_a.n_exons += 1;
        } else {
            return -1_000_008;
        }
    }

    let ne = tr_a.n_exons as usize - 1;
    tr_a.exons[ne][EX_I_FRAG] = i_frag_b;
    tr_a.exons[ne][EX_SJ_A] = sj_ab;

    // suppress unused bindings triggered by feature-gated match arms
    let _ = g_aend;
    let _ = MARK_FRAG_SPACER_BASE;

    score.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

/// Tiny helper used in the insertion case to compute `r1 = rAend + ii + (ii<=jR ? 0:Ins)`.
#[inline]
fn r_aend_plus_ins(r_aend: u64, ii: i64, j_r: i64, ins_l: i64) -> i64 {
    let pad = if ii <= j_r { 0 } else { ins_l };
    r_aend as i64 + ii + pad
}
