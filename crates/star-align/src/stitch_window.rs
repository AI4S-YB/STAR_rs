//! Port of `stitchWindowAligns.cpp` (the recursive stitcher that expands a
//! window's seed list into a set of transcripts) and
//! `ReadAlign::stitchWindowSeeds` (the driver for long reads; we keep a stub
//! that delegates to `stitch_window_aligns` under `#[cfg(feature = "long-reads")]`).

use star_core::types::{
    EX_G, EX_I_FRAG, EX_L, EX_R, EX_SJ_A, WA_ANCHOR, WA_G_START, WA_I_FRAG, WA_LENGTH, WA_NREP,
    WA_R_START, WA_SJ_A,
};
use star_core::SCORE_MATCH;
use star_genome::Genome;
use star_params::parameters::Parameters;

use crate::read_align::{ReadAlign, WinAlign};
use crate::stitch::{
    binary_search2, blocks_overlap, extend_align, stitch_align_to_transcript,
};
use crate::transcript::Transcript;

/// Snapshot of the host `ReadAlign` state that `stitch_window_aligns` is
/// allowed to mutate. The C++ code accesses `RA->maxScoreMate`,
/// `RA->outFilterMismatchNmaxTotal`, `RA->readLength` and `RA->iRead`; we
/// surface just these.
pub struct StitchCtx<'a> {
    pub max_score_mate: &'a mut [i32],
    pub out_filter_mismatch_nmax_total: u64,
    pub read_length: &'a [u64],
    pub i_read: u64,
}

/// 1:1 port of `stitchWindowAligns` (`stitchWindowAligns.cpp`).
///
/// Only single-chromosome single-fragment short-read stitching is implemented
/// in this first pass; long-reads-specific branches (COMPILE_FOR_LONG_READS)
/// are omitted since they are controlled by a build feature.
#[allow(clippy::too_many_arguments)]
pub unsafe fn stitch_window_aligns(
    i_a: u64,
    n_a: u64,
    mut score: i32,
    wa_incl: &mut [bool],
    t_r2: u64,
    t_g2: u64,
    tr_a: Transcript,
    l_read: u64,
    wa: &[WinAlign],
    r_forward: *const u8,
    map_gen: &Genome,
    p: &Parameters,
    w_tr: &mut Vec<Transcript>,
    n_win_tr: &mut u64,
    ctx: &mut StitchCtx,
) {
    let mut t_r2 = t_r2;
    let mut t_g2 = t_g2;
    let mut tr_a = tr_a;

    if i_a >= n_a && t_r2 == 0 {
        return;
    }

    if i_a >= n_a {
        // Extend (5' first if on forward read-strand, else 3' first).
        let v_order: [usize; 2] = if tr_a.ro_str == 0 { [0, 1] } else { [1, 0] };

        for &step in &v_order {
            match step {
                0 => {
                    if tr_a.r_start > 0 {
                        let mut tr_step = Transcript::new();
                        tr_step.reset();
                        let imate = tr_a.exons[0][EX_I_FRAG] as usize;
                        let end_bit = (tr_a.str_ != imate as u64) as usize;
                        if extend_align(
                            r_forward,
                            map_gen.g.as_ptr().add(star_genome::load::LOAD_L),
                            tr_a.r_start - 1,
                            tr_a.g_start - 1,
                            -1,
                            -1,
                            tr_a.r_start,
                            t_r2 - tr_a.r_start + 1,
                            tr_a.n_mm,
                            ctx.out_filter_mismatch_nmax_total,
                            p.out_filter_mismatch_nover_lmax,
                            p.align_ends_type.ext[imate][end_bit],
                            &mut tr_step,
                        ) {
                            tr_a.add(&tr_step);
                            score += tr_step.max_score;
                            tr_a.r_start -= tr_step.extend_l;
                            tr_a.g_start -= tr_step.extend_l;
                            tr_a.exons[0][EX_R] = tr_a.r_start;
                            tr_a.exons[0][EX_G] = tr_a.g_start;
                            tr_a.exons[0][EX_L] += tr_step.extend_l;
                        }
                    }
                }
                1 => {
                    if t_r2 < l_read {
                        let mut tr_step = Transcript::new();
                        tr_step.reset();
                        let last = tr_a.n_exons as usize - 1;
                        let imate = tr_a.exons[last][EX_I_FRAG] as usize;
                        let end_bit = (imate as u64 == tr_a.str_) as usize;
                        if extend_align(
                            r_forward,
                            map_gen.g.as_ptr().add(star_genome::load::LOAD_L),
                            t_r2 + 1,
                            t_g2 + 1,
                            1,
                            1,
                            l_read - t_r2 - 1,
                            t_r2 - tr_a.r_start + 1,
                            tr_a.n_mm,
                            ctx.out_filter_mismatch_nmax_total,
                            p.out_filter_mismatch_nover_lmax,
                            p.align_ends_type.ext[imate][end_bit],
                            &mut tr_step,
                        ) {
                            tr_a.add(&tr_step);
                            score += tr_step.max_score;
                            t_r2 += tr_step.extend_l;
                            t_g2 += tr_step.extend_l;
                            tr_a.exons[last][EX_L] += tr_step.extend_l;
                        }
                    }
                }
                _ => {}
            }
        }

        // Chromosome-boundary clipping
        let soft_clip = p.align_soft_clip_at_reference_ends == "Yes";
        let chr = tr_a.chr as usize;
        let last = tr_a.n_exons as usize - 1;
        if !soft_clip
            && (tr_a.exons[last][EX_G] + l_read - tr_a.exons[last][EX_R]
                > map_gen.chr_start[chr] + map_gen.chr_length[chr]
                || tr_a.exons[0][EX_G] < map_gen.chr_start[chr] + tr_a.exons[0][EX_R])
        {
            return;
        }

        tr_a.r_length = 0;
        for isj in 0..tr_a.n_exons as usize {
            tr_a.r_length += tr_a.exons[isj][EX_L];
        }
        tr_a.g_length = t_g2 + 1 - tr_a.g_start;

        // Minimum junction-overhang filter.
        for isj in 0..tr_a.n_exons.saturating_sub(1) as usize {
            if tr_a.canon_sj[isj] >= 0 {
                if tr_a.sj_annot[isj] == 1 {
                    let left_ok = tr_a.exons[isj][EX_L] >= p.align_sjdb_overhang_min
                        || !(isj == 0
                            || tr_a.canon_sj[isj - 1] == -3
                            || (tr_a.sj_annot[isj - 1] == 0 && tr_a.canon_sj[isj - 1] >= 0));
                    let right_ok = tr_a.exons[isj + 1][EX_L] >= p.align_sjdb_overhang_min
                        || !(isj + 1 == tr_a.n_exons as usize - 1
                            || tr_a.canon_sj[isj + 1] == -3
                            || (tr_a.sj_annot[isj + 1] == 0
                                && tr_a.canon_sj[isj + 1] >= 0));
                    if !(left_ok && right_ok) {
                        return;
                    }
                } else {
                    if tr_a.exons[isj][EX_L]
                        < p.align_sj_overhang_min + tr_a.shift_sj[isj][0]
                        || tr_a.exons[isj + 1][EX_L]
                            < p.align_sj_overhang_min + tr_a.shift_sj[isj][1]
                    {
                        return;
                    }
                }
            }
        }
        if tr_a.n_exons > 1
            && tr_a.sj_annot[tr_a.n_exons as usize - 2] == 1
            && tr_a.exons[tr_a.n_exons as usize - 1][EX_L] < p.align_sjdb_overhang_min
        {
            return;
        }

        // Strand consistency
        let mut sj_n: u64 = 0;
        tr_a.intron_motifs = [0; 3];
        tr_a.sj_yes = false;
        for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
            if tr_a.canon_sj[iex] >= 0 {
                sj_n += 1;
                let k = tr_a.sj_str[iex] as usize;
                if k < 3 {
                    tr_a.intron_motifs[k] += 1;
                }
                tr_a.sj_yes = true;
            }
        }
        tr_a.sj_motif_strand = if tr_a.intron_motifs[1] > 0 && tr_a.intron_motifs[2] == 0 {
            1
        } else if tr_a.intron_motifs[1] == 0 && tr_a.intron_motifs[2] > 0 {
            2
        } else {
            0
        };
        if tr_a.intron_motifs[1] > 0
            && tr_a.intron_motifs[2] > 0
            && p.out_filter_intron_strands == "RemoveInconsistentStrands"
        {
            return;
        }
        if sj_n > 0 && tr_a.sj_motif_strand == 0 && p.out_sam_strand_field_type == 1 {
            return;
        }

        match p.out_filter_intron_motifs.as_str() {
            "None" => {}
            "RemoveNoncanonical" => {
                for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
                    if tr_a.canon_sj[iex] == 0 {
                        return;
                    }
                }
            }
            "RemoveNoncanonicalUnannotated" => {
                for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
                    if tr_a.canon_sj[iex] == 0 && tr_a.sj_annot[iex] == 0 {
                        return;
                    }
                }
            }
            _ => {} // defer error to finalize
        }

        // Mapped-length-per-mate check
        {
            let mut nsj: u64 = 0;
            let mut exl: u64 = 0;
            for iex in 0..tr_a.n_exons as usize {
                exl += tr_a.exons[iex][EX_L];
                let last_or_mate_end = iex == tr_a.n_exons as usize - 1
                    || tr_a.canon_sj[iex] == -3;
                if last_or_mate_end {
                    // C++ stitchWindowAligns.cpp:159 casts
                    // `(uint)(alignSplicedMateMapLminOverLmate * readLength[...])`
                    // before comparing `exl < ...`. Preserve the truncation so
                    // that e.g. `exl=99` passes when the threshold is 0.66*151.
                    let lmate = ctx.read_length[tr_a.exons[iex][EX_I_FRAG] as usize];
                    let lmate_frac =
                        (p.align_spliced_mate_map_lmin_over_lmate * lmate as f64) as u64;
                    if nsj > 0
                        && (exl < p.align_spliced_mate_map_lmin || exl < lmate_frac)
                    {
                        return;
                    }
                    exl = 0;
                    nsj = 0;
                } else if tr_a.canon_sj[iex] >= 0 {
                    nsj += 1;
                }
            }
        }

        if p.out_filter_by_sjout_stage == 2 {
            for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
                if tr_a.canon_sj[iex] >= 0 && tr_a.sj_annot[iex] == 0 {
                    let j_s = tr_a.exons[iex][EX_G] + tr_a.exons[iex][EX_L];
                    let j_e = tr_a.exons[iex + 1][EX_G] - 1;
                    if binary_search2(j_s, j_e, &p.sj_novel_start, &p.sj_novel_end) < 0 {
                        return;
                    }
                }
            }
        }

        // Mate-overlap junction consistency (paired-end only — single-end
        // paths will early-exit this block because iFrags match).
        if tr_a.exons[0][EX_I_FRAG] != tr_a.exons[tr_a.n_exons as usize - 1][EX_I_FRAG] {
            let last = tr_a.n_exons as usize - 1;
            if tr_a.exons[last][EX_G] + tr_a.exons[last][EX_L] <= tr_a.exons[0][EX_G] {
                return;
            }
            let mut iex_m2 = tr_a.n_exons as usize;
            for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
                if tr_a.canon_sj[iex] == -3 {
                    iex_m2 = iex + 1;
                    break;
                }
            }
            if iex_m2 < tr_a.n_exons as usize
                && tr_a.exons[iex_m2 - 1][EX_G] + tr_a.exons[iex_m2 - 1][EX_L]
                    > tr_a.exons[iex_m2][EX_G]
            {
                if tr_a.exons[0][EX_G]
                    > tr_a.exons[iex_m2][EX_G]
                        + tr_a.exons[0][EX_R]
                        + p.align_ends_protrude.n_bases_max as u64
                {
                    return;
                }
                if tr_a.exons[iex_m2 - 1][EX_G] + tr_a.exons[iex_m2 - 1][EX_L]
                    > tr_a.exons[last][EX_G] + l_read - tr_a.exons[last][EX_R]
                        + p.align_ends_protrude.n_bases_max as u64
                {
                    return;
                }
                let mut iex1 = 1;
                let mut iex2 = iex_m2 + 1;
                while iex1 < iex_m2 {
                    if tr_a.exons[iex1][EX_G]
                        >= tr_a.exons[iex2 - 1][EX_G] + tr_a.exons[iex2 - 1][EX_L]
                    {
                        break;
                    }
                    iex1 += 1;
                }
                while iex1 < iex_m2 && iex2 < tr_a.n_exons as usize {
                    if tr_a.canon_sj[iex1 - 1] < 0 {
                        iex1 += 1;
                        continue;
                    }
                    if tr_a.canon_sj[iex2 - 1] < 0 {
                        iex2 += 1;
                        continue;
                    }
                    if tr_a.exons[iex1][EX_G] != tr_a.exons[iex2][EX_G]
                        || (tr_a.exons[iex1 - 1][EX_G] + tr_a.exons[iex1 - 1][EX_L])
                            != (tr_a.exons[iex2 - 1][EX_G] + tr_a.exons[iex2 - 1][EX_L])
                    {
                        return;
                    }
                    iex1 += 1;
                    iex2 += 1;
                }
            }
        }

        if p.score_genomic_length_log2_scale != 0.0 {
            let last = tr_a.n_exons as usize - 1;
            let glen = (tr_a.exons[last][EX_G] + tr_a.exons[last][EX_L] - tr_a.exons[0][EX_G])
                as f64;
            score += (glen.log2() * p.score_genomic_length_log2_scale - 0.5).ceil() as i32;
            if score < 0 {
                score = 0;
            }
        }

        tr_a.ro_start = if tr_a.ro_str == 0 {
            tr_a.r_start
        } else {
            l_read - tr_a.r_start - tr_a.r_length
        };
        tr_a.max_score = score;

        if tr_a.exons[0][EX_I_FRAG] == tr_a.exons[tr_a.n_exons as usize - 1][EX_I_FRAG] {
            tr_a.i_frag = tr_a.exons[0][EX_I_FRAG] as i32;
            let f = tr_a.i_frag as usize;
            if f < ctx.max_score_mate.len() && score > ctx.max_score_mate[f] {
                ctx.max_score_mate[f] = score;
            }
        } else {
            tr_a.i_frag = -1;
        }

        // `variationAdjust` — not yet ported; returns 0 here.
        tr_a.max_score = score;

        let max_score_gate = if w_tr.is_empty() {
            0
        } else {
            w_tr[0].max_score
        };
        let pass = score + p.out_filter_multimap_score_range as i32 >= max_score_gate
            || (tr_a.i_frag >= 0
                && score + p.out_filter_multimap_score_range as i32
                    >= ctx.max_score_mate[tr_a.i_frag as usize])
            || p.p_ch_segment_min > 0;


        if pass {
            tr_a.mapped_length = tr_a.exons.iter().take(tr_a.n_exons as usize).map(|e| e[EX_L]).sum();

            let mut i_tr = 0usize;
            while i_tr < *n_win_tr as usize {
                let n_overlap = blocks_overlap(&tr_a, &w_tr[i_tr]);
                let u_new = tr_a.mapped_length.saturating_sub(n_overlap);
                let u_old = w_tr[i_tr].mapped_length.saturating_sub(n_overlap);
                if u_new == 0 && score < w_tr[i_tr].max_score {
                    break;
                } else if u_old == 0 {
                    w_tr.remove(i_tr);
                    *n_win_tr -= 1;
                } else if u_old > 0 && (u_new > 0 || score >= w_tr[i_tr].max_score) {
                    i_tr += 1;
                } else {
                    break;
                }
            }
            if i_tr == *n_win_tr as usize {
                let mut pos = 0usize;
                while pos < *n_win_tr as usize {
                    if score > w_tr[pos].max_score
                        || (score == w_tr[pos].max_score
                            && tr_a.g_length < w_tr[pos].g_length)
                    {
                        break;
                    }
                    pos += 1;
                }
                w_tr.insert(pos, tr_a.clone());
                if *n_win_tr < p.align_transcripts_per_window_nmax as u64 {
                    *n_win_tr += 1;
                }
            }
        }
        return;
    }

    // Recursive branch — include or exclude current align.
    let d_score: i32;
    let mut tr_ai = tr_a.clone();
    if tr_a.n_exons > 0 {
        d_score = stitch_align_to_transcript(
            t_r2,
            t_g2,
            wa[i_a as usize][WA_R_START],
            wa[i_a as usize][WA_G_START],
            wa[i_a as usize][WA_LENGTH],
            wa[i_a as usize][WA_I_FRAG],
            wa[i_a as usize][WA_SJ_A],
            p,
            r_forward,
            map_gen,
            &mut tr_ai,
            ctx.out_filter_mismatch_nmax_total,
        );
    } else {
        tr_ai.exons[0][EX_R] = wa[i_a as usize][WA_R_START];
        tr_ai.r_start = wa[i_a as usize][WA_R_START];
        tr_ai.exons[0][EX_G] = wa[i_a as usize][WA_G_START];
        tr_ai.g_start = wa[i_a as usize][WA_G_START];
        tr_ai.exons[0][EX_L] = wa[i_a as usize][WA_LENGTH];
        tr_ai.exons[0][EX_I_FRAG] = wa[i_a as usize][WA_I_FRAG];
        tr_ai.exons[0][EX_SJ_A] = wa[i_a as usize][WA_SJ_A];
        tr_ai.n_exons = 1;
        d_score = (wa[i_a as usize][WA_LENGTH] as i64 * SCORE_MATCH as i64) as i32;
        tr_ai.n_match = wa[i_a as usize][WA_LENGTH];
        for e in wa_incl.iter_mut() {
            *e = false;
        }
    }

    if d_score > -1_000_000 {
        wa_incl[i_a as usize] = true;
        if wa[i_a as usize][WA_NREP] == 1 {
            tr_ai.n_unique += 1;
        }
        if wa[i_a as usize][WA_ANCHOR] > 0 {
            tr_ai.n_anchor += 1;
        }
        stitch_window_aligns(
            i_a + 1,
            n_a,
            score + d_score,
            wa_incl,
            wa[i_a as usize][WA_R_START] + wa[i_a as usize][WA_LENGTH] - 1,
            wa[i_a as usize][WA_G_START] + wa[i_a as usize][WA_LENGTH] - 1,
            tr_ai,
            l_read,
            wa,
            r_forward,
            map_gen,
            p,
            w_tr,
            n_win_tr,
            ctx,
        );
    }

    if wa[i_a as usize][WA_ANCHOR] != 2 || tr_a.n_anchor > 0 {
        wa_incl[i_a as usize] = false;
        stitch_window_aligns(
            i_a + 1,
            n_a,
            score,
            wa_incl,
            t_r2,
            t_g2,
            tr_a,
            l_read,
            wa,
            r_forward,
            map_gen,
            p,
            w_tr,
            n_win_tr,
            ctx,
        );
    }
}

impl ReadAlign {
    /// Stub for `stitchWindowSeeds` — only used in COMPILE_FOR_LONG_READS
    /// builds. Delegates to `stitch_window_aligns` in short-read mode.
    pub unsafe fn stitch_window_seeds(
        &mut self,
        _p: &Parameters,
        _map_gen: &Genome,
        _i_w: u64,
        _i_w_rec: u64,
        _wa_excl: Option<&[bool]>,
        _r_forward: *const u8,
    ) {
        // Short-read path uses `stitch_window_aligns` directly from
        // `stitch_pieces`; long-read path to be ported in a later milestone.
    }
}
