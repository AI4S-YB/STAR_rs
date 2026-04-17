//! 1:1 port of `ReadAlign::chimericDetection` and
//! `ReadAlign::chimericDetectionOld` (ReadAlign_chimericDetection.cpp and
//! ReadAlign_chimericDetectionOld.cpp).
//!
//! Only the "Old" legacy path is implemented here. The "Mult" path
//! (`ChimericDetection::chimericDetectionMult`, enabled by
//! `--chimMultimapNmax>0`) and the PE-merged path
//! (`ReadAlign::chimericDetectionPEmerged`) are stubbed and fall through to
//! "no chimeric detected" until they are wired up in a later sub-milestone.

use std::io::Write;

use star_align::read_align::ReadAlign;
use star_align::stitch::blocks_overlap;
use star_align::transcript::Transcript;
use star_core::{EX_G, EX_I_FRAG, EX_L, EX_R};
use star_genome::Genome;
use star_genome::load::LOAD_L;
use star_params::parameters::Parameters;

/// Port of `ReadAlign::chimericDetection` (entry point).
///
/// Dispatches to the Old detector when `--chimMultimapNmax==0`, and
/// records `ra.chim_record` on success. Returns `Ok(())` when no error.
pub fn chimeric_detection(
    ra: &mut ReadAlign,
    p: &Parameters,
    map_gen: &Genome,
) -> anyhow::Result<()> {
    ra.chim_record = false;
    if p.p_ch.segment_min == 0 {
        return Ok(());
    }
    if p.out_filter_by_sjout_stage > 1 {
        // No chimeric output for stage=2.
        return Ok(());
    }

    if p.p_ch.multimap_nmax == 0 {
        ra.chim_record = chimeric_detection_old(ra, p, map_gen);
    } else {
        // chimericDetectionMult path — not implemented yet.
        // C++: ChimericDetection::chimericDetectionMult(nW, readLength,
        //      trBest->maxScore, NULL).
        ra.chim_record = false;
    }

    Ok(())
}

/// 1:1 port of `ReadAlign::chimericDetectionOld()`.
///
/// Returns true if a chimeric alignment was recorded (trChim[0..1],
/// chimJ0/J1, chimMotif, chimStr, chimRepeat0/1).
pub fn chimeric_detection_old(
    ra: &mut ReadAlign,
    p: &Parameters,
    map_gen: &Genome,
) -> bool {
    let p_ch = &p.p_ch;
    let read_length0 = ra.read_length[0];
    let read_length1 = ra.read_length.get(1).copied().unwrap_or(0);
    let l_read = ra.l_read;

    if ra.n_tr > p_ch.main_segment_mult_nmax && ra.n_tr != 2 {
        return false;
    }

    let tr_best_idx = ra.tr_best as usize;
    if tr_best_idx >= ra.tr_all.len() || ra.tr_all[tr_best_idx].is_empty() {
        return false;
    }

    // Snapshot trBest as a local struct (we need to mutate ra.tr_all below
    // after trChim[0]=*trBest).
    let tr_best_snapshot: Transcript = ra.tr_all[tr_best_idx][0].clone();

    let last_best = tr_best_snapshot.n_exons as usize - 1;
    let tr_best_cov_end = tr_best_snapshot.exons[last_best][EX_R]
        + tr_best_snapshot.exons[last_best][EX_L];
    if !(p_ch.segment_min > 0
        && tr_best_snapshot.r_length >= p_ch.segment_min
        && (tr_best_cov_end + p_ch.segment_min <= l_read
            || tr_best_snapshot.exons[0][EX_R] >= p_ch.segment_min)
        && tr_best_snapshot.intron_motifs[0] == 0
        && (tr_best_snapshot.intron_motifs[1] == 0 || tr_best_snapshot.intron_motifs[2] == 0))
    {
        return false;
    }

    let mut chim_score_best: i32 = 0;
    let mut chim_score_next: i32 = 0;
    ra.tr_chim[0] = tr_best_snapshot.clone();

    let mut ro_start1 = if tr_best_snapshot.str_ == 0 {
        tr_best_snapshot.exons[0][EX_R]
    } else {
        l_read
            - tr_best_snapshot.exons[last_best][EX_R]
            - tr_best_snapshot.exons[last_best][EX_L]
    };
    let mut ro_end1 = if tr_best_snapshot.str_ == 0 {
        tr_best_snapshot.exons[last_best][EX_R] + tr_best_snapshot.exons[last_best][EX_L] - 1
    } else {
        l_read - tr_best_snapshot.exons[0][EX_R] - 1
    };
    if ro_start1 > read_length0 {
        ro_start1 -= 1;
    }
    if ro_end1 > read_length0 {
        ro_end1 -= 1;
    }

    let mut chim_str_best: u64 = 0;
    let mut chim_str: u64 = if tr_best_snapshot.intron_motifs[1] == 0
        && tr_best_snapshot.intron_motifs[2] == 0
    {
        0
    } else if (tr_best_snapshot.str_ == 0) == (tr_best_snapshot.intron_motifs[1] > 0) {
        1
    } else {
        2
    };

    // Index of the chosen trChim[1] in (iW, iWt) — used for the mainSegmentMultNmax==2 check.
    let mut tr_chim1_idx: Option<(usize, usize)> = None;
    // tr_chim[1] pre-assignment snapshot, used for blocksOverlap between candidates.
    let mut tr_chim1_prev: Option<Transcript> = None;

    let n_w = ra.n_w as usize;
    for i_w in 0..n_w {
        let n_win_tr = ra.n_win_tr[i_w] as usize;
        for i_wt in 0..n_win_tr {
            // Identify whether this window's best is the global best.
            let i_w_best = {
                // `trBest != trAll[iW][0]` — we compare by cloning pointer
                // semantics: the original code uses pointer equality, which
                // is true iff this is the best window (tr_best index == i_w).
                i_w == tr_best_idx
            };
            if !i_w_best && i_wt > 0 {
                break;
            }
            if i_w_best && i_wt == 0 {
                continue;
            }
            let cand = &ra.tr_all[i_w][i_wt];
            if cand.intron_motifs[0] > 0 {
                continue;
            }
            let chim_str1: u64 = if cand.intron_motifs[1] == 0 && cand.intron_motifs[2] == 0 {
                0
            } else if (cand.str_ == 0) == (cand.intron_motifs[1] > 0) {
                1
            } else {
                2
            };
            if chim_str != 0 && chim_str1 != 0 && chim_str != chim_str1 {
                continue;
            }

            let last_cand = cand.n_exons as usize - 1;
            let mut ro_start2 = if cand.str_ == 0 {
                cand.exons[0][EX_R]
            } else {
                l_read - cand.exons[last_cand][EX_R] - cand.exons[last_cand][EX_L]
            };
            let mut ro_end2 = if cand.str_ == 0 {
                cand.exons[last_cand][EX_R] + cand.exons[last_cand][EX_L] - 1
            } else {
                l_read - cand.exons[0][EX_R] - 1
            };
            if ro_start2 > read_length0 {
                ro_start2 -= 1;
            }
            if ro_end2 > read_length0 {
                ro_end2 -= 1;
            }

            let chim_overlap: u64 = if ro_start2 > ro_start1 {
                if ro_start2 > ro_end1 {
                    0
                } else {
                    ro_end1 - ro_start2 + 1
                }
            } else if ro_end2 < ro_start1 {
                0
            } else {
                ro_end2 - ro_start1 + 1
            };
            let diff_mates = (ro_end1 < read_length0 && ro_start2 >= read_length0)
                || (ro_end2 < read_length0 && ro_start1 >= read_length0);

            if ro_end1 > p_ch.segment_min + ro_start1 + chim_overlap
                && ro_end2 > p_ch.segment_min + ro_start2 + chim_overlap
                && (diff_mates
                    || ((ro_end1 + p_ch.segment_read_gap_max + 1) >= ro_start2
                        && (ro_end2 + p_ch.segment_read_gap_max + 1) >= ro_start1))
            {
                let chim_score: i32 = tr_best_snapshot.max_score + cand.max_score
                    - chim_overlap as i32;

                let overlap1 = if i_wt > 0 && chim_score_best > 0 {
                    if let Some(prev) = &tr_chim1_prev {
                        blocks_overlap(prev, cand)
                    } else {
                        0
                    }
                } else {
                    0
                };

                if chim_score > chim_score_best {
                    // Record previous tr_chim[1] for the overlap check on the
                    // next iteration; match pointer semantics in C++ where
                    // `trChim[1]` holds the previous candidate before being
                    // overwritten.
                    tr_chim1_prev = Some(cand.clone());
                    ra.tr_chim[1] = cand.clone();
                    tr_chim1_idx = Some((i_w, i_wt));
                    if overlap1 == 0 {
                        chim_score_next = chim_score_best;
                    }
                    chim_score_best = chim_score;
                    // Update roStart / cStart on trChim[1] (C++ lines 89-90).
                    let tc1 = &mut ra.tr_chim[1];
                    tc1.ro_start = if tc1.ro_str == 0 {
                        tc1.r_start
                    } else {
                        l_read - tc1.r_start - tc1.r_length
                    };
                    tc1.c_start = tc1
                        .g_start
                        .saturating_sub(map_gen.chr_start[tc1.chr as usize]);
                    chim_str_best = chim_str1;
                } else if chim_score > chim_score_next && overlap1 == 0 {
                    chim_score_next = chim_score;
                }
            }
        }
    }

    if !(chim_score_best >= p_ch.score_min
        && chim_score_best + p_ch.score_drop_max >= (read_length0 + read_length1) as i32)
    {
        return false;
    }

    // nTr>mainSegmentMultNmax case (nTr==2 special).
    if ra.n_tr > p_ch.main_segment_mult_nmax {
        if let Some((iw, iwt)) = tr_chim1_idx {
            // C++ compares trChim1 pointer against trMult[0] / trMult[1]. In
            // our port trMult contains window-indices; check that the
            // chosen (iw,iwt) corresponds to one of the first two multi
            // transcripts. Since multMapSelect stores transcripts in
            // tr_mult_array by order of (iw,iwt), we compare by content.
            let pick = &ra.tr_all[iw][iwt];
            let mut ok = false;
            for k in 0..ra.tr_mult_array.len().min(2) {
                let tm = &ra.tr_mult_array[k];
                if tm.g_start == pick.g_start
                    && tm.r_start == pick.r_start
                    && tm.str_ == pick.str_
                    && tm.chr == pick.chr
                    && tm.n_exons == pick.n_exons
                {
                    ok = true;
                    break;
                }
            }
            if !ok {
                return false;
            }
        } else {
            return false;
        }
    }

    if chim_str == 0 {
        chim_str = chim_str_best;
    }

    ra.chim_n = 0;
    if chim_score_next + p_ch.score_separation >= chim_score_best {
        return false;
    }
    if ra.tr_chim[0].ro_start > ra.tr_chim[1].ro_start {
        ra.tr_chim.swap(0, 1);
    }

    let e0 = if ra.tr_chim[0].str_ == 1 {
        0usize
    } else {
        ra.tr_chim[0].n_exons as usize - 1
    };
    let e1 = if ra.tr_chim[1].str_ == 0 {
        0usize
    } else {
        ra.tr_chim[1].n_exons as usize - 1
    };

    let mut chim_repeat0: u64 = 0;
    let mut chim_repeat1: u64 = 0;
    let chim_j0: u64;
    let chim_j1: u64;
    let mut chim_motif: i32 = 0;

    ra.chim_n = 2;
    let tr0_i_frag = ra.tr_chim[0].exons[e0][EX_I_FRAG];
    let tr1_i_frag = ra.tr_chim[1].exons[e1][EX_I_FRAG];

    if tr0_i_frag > tr1_i_frag {
        // Strange configuration — rejected.
        return false;
    } else if tr0_i_frag < tr1_i_frag {
        // Mates bracket the chimeric junction.
        ra.chim_n = 2;
        chim_motif = -1;
        if ra.tr_chim[0].str_ == 1 {
            chim_j0 = ra.tr_chim[0].exons[e0][EX_G] - 1;
        } else {
            chim_j0 = ra.tr_chim[0].exons[e0][EX_G] + ra.tr_chim[0].exons[e0][EX_L];
        }
        if ra.tr_chim[1].str_ == 0 {
            chim_j1 = ra.tr_chim[1].exons[e1][EX_G] - 1;
        } else {
            chim_j1 = ra.tr_chim[1].exons[e1][EX_G] + ra.tr_chim[1].exons[e1][EX_L];
        }
    } else {
        // Chimeric junction is within one of the mates; need to scan for motif.
        if !(ra.tr_chim[0].exons[e0][EX_L] >= p_ch.junction_overhang_min
            && ra.tr_chim[1].exons[e1][EX_L] >= p_ch.junction_overhang_min)
        {
            return false;
        }
        let ro_start0 = if ra.tr_chim[0].str_ == 0 {
            ra.tr_chim[0].exons[e0][EX_R]
        } else {
            l_read - ra.tr_chim[0].exons[e0][EX_R] - ra.tr_chim[0].exons[e0][EX_L]
        };
        let ro_start1_in = if ra.tr_chim[1].str_ == 0 {
            ra.tr_chim[1].exons[e1][EX_R]
        } else {
            l_read - ra.tr_chim[1].exons[e1][EX_R] - ra.tr_chim[1].exons[e1][EX_L]
        };

        let mut j_r_best: u64 = 0;
        let mut j_score: i32 = 0;
        let mut j_score_best: i32 = -999_999;

        let mut j_r_max = ro_start1_in + ra.tr_chim[1].exons[e1][EX_L];
        j_r_max = if j_r_max > ro_start0 {
            j_r_max - ro_start0 - 1
        } else {
            0
        };

        let read1_fwd = &ra.read1[0];
        let mut j_r: u64 = 0;
        while j_r < j_r_max {
            if j_r == read_length0 {
                j_r += 1;
            }
            if j_r >= j_r_max {
                break;
            }
            let b_r = read1_fwd[(ro_start0 + j_r) as usize];

            let b0 = if ra.tr_chim[0].str_ == 0 {
                map_gen.g[(ra.tr_chim[0].exons[e0][EX_G] + j_r) as usize + LOAD_L]
            } else {
                let raw = map_gen.g[(ra.tr_chim[0].exons[e0][EX_G]
                    + ra.tr_chim[0].exons[e0][EX_L]
                    - 1
                    - j_r) as usize + LOAD_L];
                if raw < 4 { 3 - raw } else { raw }
            };
            let b1 = if ra.tr_chim[1].str_ == 0 {
                map_gen.g[(ra.tr_chim[1].exons[e1][EX_G] - ro_start1_in + ro_start0 + j_r)
                    as usize + LOAD_L]
            } else {
                let raw = map_gen.g[(ra.tr_chim[1].exons[e1][EX_G]
                    + ra.tr_chim[1].exons[e1][EX_L]
                    - 1
                    + ro_start1_in
                    - ro_start0
                    - j_r) as usize + LOAD_L];
                if raw < 4 { 3 - raw } else { raw }
            };

            if (p_ch.filter_genomic_n && (b0 > 3 || b1 > 3)) || b_r > 3 {
                ra.chim_n = 0;
                break;
            }

            let (b01, b02) = if ra.tr_chim[0].str_ == 0 {
                (
                    map_gen.g[(ra.tr_chim[0].exons[e0][EX_G] + j_r + 1) as usize + LOAD_L],
                    map_gen.g[(ra.tr_chim[0].exons[e0][EX_G] + j_r + 2) as usize + LOAD_L],
                )
            } else {
                let r01 = map_gen.g[(ra.tr_chim[0].exons[e0][EX_G]
                    + ra.tr_chim[0].exons[e0][EX_L]
                    - 1
                    - j_r
                    - 1) as usize + LOAD_L];
                let r02 = map_gen.g[(ra.tr_chim[0].exons[e0][EX_G]
                    + ra.tr_chim[0].exons[e0][EX_L]
                    - 1
                    - j_r
                    - 2) as usize + LOAD_L];
                (
                    if r01 < 4 { 3 - r01 } else { r01 },
                    if r02 < 4 { 3 - r02 } else { r02 },
                )
            };
            let (b11, b12) = if ra.tr_chim[1].str_ == 0 {
                (
                    map_gen.g[(ra.tr_chim[1].exons[e1][EX_G] - ro_start1_in + ro_start0 + j_r
                        - 1) as usize + LOAD_L],
                    map_gen.g[(ra.tr_chim[1].exons[e1][EX_G] - ro_start1_in + ro_start0 + j_r)
                        as usize + LOAD_L],
                )
            } else {
                let r11 = map_gen.g[(ra.tr_chim[1].exons[e1][EX_G]
                    + ra.tr_chim[1].exons[e1][EX_L]
                    - 1
                    + ro_start1_in
                    - ro_start0
                    - j_r
                    + 1) as usize + LOAD_L];
                let r12 = map_gen.g[(ra.tr_chim[1].exons[e1][EX_G]
                    + ra.tr_chim[1].exons[e1][EX_L]
                    - 1
                    + ro_start1_in
                    - ro_start0
                    - j_r) as usize + LOAD_L];
                (
                    if r11 < 4 { 3 - r11 } else { r11 },
                    if r12 < 4 { 3 - r12 } else { r12 },
                )
            };

            let mut j_motif: i32 = 0;
            if b01 == 2 && b02 == 3 && b11 == 0 && b12 == 2 {
                if chim_str != 2 {
                    j_motif = 1;
                }
            } else if b01 == 1 && b02 == 3 && b11 == 0 && b12 == 1 {
                if chim_str != 1 {
                    j_motif = 2;
                }
            }

            if b_r == b0 && b_r != b1 {
                j_score += 1;
            } else if b_r != b0 && b_r == b1 {
                j_score -= 1;
            }

            let j_score_j = if j_motif == 0 {
                j_score + p_ch.score_junction_non_gtag
            } else {
                j_score
            };
            if j_score_j > j_score_best || (j_score_j == j_score_best && j_motif > 0) {
                chim_motif = j_motif;
                j_r_best = j_r;
                j_score_best = j_score_j;
            }

            j_r += 1;
        }
        if ra.chim_n == 0 {
            return false;
        }

        if chim_motif == 0 {
            let adj = chim_score_best + 1 + p_ch.score_junction_non_gtag;
            chim_score_best = adj;
            if !(chim_score_best >= p_ch.score_min
                && chim_score_best + p_ch.score_drop_max
                    >= (read_length0 + read_length1) as i32)
            {
                return false;
            }
        }

        // Shift junction in trChim.
        if ra.tr_chim[0].str_ == 1 {
            let inc = ra.tr_chim[0].exons[e0][EX_L] - j_r_best - 1;
            ra.tr_chim[0].exons[e0][EX_R] += inc;
            ra.tr_chim[0].exons[e0][EX_G] += inc;
            ra.tr_chim[0].exons[e0][EX_L] = j_r_best + 1;
            chim_j0 = ra.tr_chim[0].exons[e0][EX_G] - 1;
        } else {
            ra.tr_chim[0].exons[e0][EX_L] = j_r_best + 1;
            chim_j0 = ra.tr_chim[0].exons[e0][EX_G] + ra.tr_chim[0].exons[e0][EX_L];
        }

        if ra.tr_chim[1].str_ == 0 {
            let delta = ro_start0 + j_r_best + 1 - ro_start1_in;
            ra.tr_chim[1].exons[e1][EX_R] += delta;
            ra.tr_chim[1].exons[e1][EX_G] += delta;
            ra.tr_chim[1].exons[e1][EX_L] =
                ro_start1_in + ra.tr_chim[1].exons[e1][EX_L] - ro_start0 - j_r_best - 1;
            chim_j1 = ra.tr_chim[1].exons[e1][EX_G] - 1;
        } else {
            ra.tr_chim[1].exons[e1][EX_L] =
                ro_start1_in + ra.tr_chim[1].exons[e1][EX_L] - ro_start0 - j_r_best - 1;
            chim_j1 = ra.tr_chim[1].exons[e1][EX_G] + ra.tr_chim[1].exons[e1][EX_L];
        }

        // Find repeats: forward.
        let mut j_r: u64 = 0;
        while j_r < 100 {
            let b0 = if ra.tr_chim[0].str_ == 0 {
                map_gen.g[(chim_j0 + j_r) as usize + LOAD_L]
            } else {
                let raw = map_gen.g[(chim_j0 - j_r) as usize + LOAD_L];
                if raw < 4 { 3 - raw } else { raw }
            };
            let b1 = if ra.tr_chim[1].str_ == 0 {
                map_gen.g[(chim_j1 + 1 + j_r) as usize + LOAD_L]
            } else {
                let raw = map_gen.g[(chim_j1 - 1 - j_r) as usize + LOAD_L];
                if raw < 4 { 3 - raw } else { raw }
            };
            if b0 != b1 {
                break;
            }
            j_r += 1;
        }
        chim_repeat1 = j_r;
        // Reverse.
        let mut j_r: u64 = 0;
        while j_r < 100 {
            let b0 = if ra.tr_chim[0].str_ == 0 {
                map_gen.g[(chim_j0 - 1 - j_r) as usize + LOAD_L]
            } else {
                let raw = map_gen.g[(chim_j0 + 1 + j_r) as usize + LOAD_L];
                if raw < 4 { 3 - raw } else { raw }
            };
            let b1 = if ra.tr_chim[1].str_ == 0 {
                map_gen.g[(chim_j1 - j_r) as usize + LOAD_L]
            } else {
                let raw = map_gen.g[(chim_j1 + j_r) as usize + LOAD_L];
                if raw < 4 { 3 - raw } else { raw }
            };
            if b0 != b1 {
                break;
            }
            j_r += 1;
        }
        chim_repeat0 = j_r;
    }

    // Final inter-chromosome / intra-chromosome gate.
    let same_chr_str = ra.tr_chim[0].str_ == ra.tr_chim[1].str_
        && ra.tr_chim[0].chr == ra.tr_chim[1].chr;
    let gap: u64 = if same_chr_str {
        if ra.tr_chim[0].str_ == 0 {
            chim_j1.saturating_sub(chim_j0) + 1
        } else {
            chim_j0.saturating_sub(chim_j1) + 1
        }
    } else {
        u64::MAX
    };
    let limit = if chim_motif >= 0 {
        p.align_intron_max
    } else {
        p.align_mates_gap_max
    };
    let far_enough = !same_chr_str || gap > limit;

    if far_enough {
        if chim_motif >= 0
            && (ra.tr_chim[0].exons[e0][EX_L] < p_ch.junction_overhang_min + chim_repeat0
                || ra.tr_chim[1].exons[e1][EX_L] < p_ch.junction_overhang_min + chim_repeat1)
        {
            return false;
        }
        ra.chim_j0 = chim_j0;
        ra.chim_j1 = chim_j1;
        ra.chim_motif = chim_motif;
        ra.chim_str = chim_str;
        ra.chim_repeat0 = chim_repeat0;
        ra.chim_repeat1 = chim_repeat1;
        return true;
    }

    false
}

/// 1:1 port of `ReadAlign::chimericDetectionOldOutput`
/// (ReadAlign_chimericDetectionOldOutput.cpp).
///
/// Writes one line to the chimeric junction stream (if `--chimOutType
/// Junctions` is active). The optional `SeparateSAMold` SAM output is not
/// yet implemented (M6 scope focuses on the junctions file).
pub fn chimeric_detection_old_output(
    ra: &mut ReadAlign,
    p: &Parameters,
    map_gen: &Genome,
    chim_junction_out: &mut dyn Write,
) -> anyhow::Result<()> {
    if !ra.chim_record {
        return Ok(());
    }

    // Re-score both chimeric pieces.
    for k in 0..2 {
        ra.tr_chim[k].align_score(
            &ra.read1,
            &map_gen.g,
            p.p_ge.sjdb_score,
            p.score_ins_base,
            p.score_ins_open,
            p.score_del_base,
            p.score_del_open,
            p.score_gap,
            p.score_gap_gcag,
            p.score_gap_atac,
            p.score_gap_noncan,
            p.score_genomic_length_log2_scale,
        );
    }

    if p.p_ch.out_junctions {
        let chr0 = &map_gen.chr_name[ra.tr_chim[0].chr as usize];
        let chr1 = &map_gen.chr_name[ra.tr_chim[1].chr as usize];
        let j0_1based = ra.chim_j0 - map_gen.chr_start[ra.tr_chim[0].chr as usize] + 1;
        let j1_1based = ra.chim_j1 - map_gen.chr_start[ra.tr_chim[1].chr as usize] + 1;
        let s0 = if ra.tr_chim[0].str_ == 0 { "+" } else { "-" };
        let s1 = if ra.tr_chim[1].str_ == 0 { "+" } else { "-" };
        let read_name = ra
            .read_name
            .strip_prefix('@')
            .or_else(|| ra.read_name.strip_prefix('>'))
            .unwrap_or(&ra.read_name);
        let start0 = ra.tr_chim[0].exons[0][EX_G] + 1
            - map_gen.chr_start[ra.tr_chim[0].chr as usize];
        let start1 = ra.tr_chim[1].exons[0][EX_G] + 1
            - map_gen.chr_start[ra.tr_chim[1].chr as usize];
        // Clone to avoid borrow conflict with tr_chim.
        let tr0 = ra.tr_chim[0].clone();
        let tr1 = ra.tr_chim[1].clone();
        let cigar0 = ra.output_transcript_cigar_p(&tr0, p);
        let cigar1 = ra.output_transcript_cigar_p(&tr1, p);

        write!(
            chim_junction_out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chr0,
            j0_1based,
            s0,
            chr1,
            j1_1based,
            s1,
            ra.chim_motif,
            ra.chim_repeat0,
            ra.chim_repeat1,
            read_name,
            start0,
            cigar0,
            start1,
            cigar1,
        )?;

        // `outSAMattrRG` is not yet ported (M6 scope). Emitted as a simple
        // newline terminator — matches STAR when `--outSAMattrRGline -`.
        writeln!(chim_junction_out)?;
    }
    Ok(())
}
