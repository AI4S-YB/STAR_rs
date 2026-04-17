//! 1:1 port of `outputSJ.cpp` — merge per-chunk `OutSJ` collectors into a
//! global list, apply SJ-filtering, and write `SJ.out.tab`.
//!
//! Filtering rules mirror `Parameters::outSJfilter*`:
//! - annotated junctions are always kept
//! - otherwise count_unique >= filter_min OR (count_unique + count_multiple)
//!   >= filter_count_total AND both overhangs >= filter_overhang_min.
//! Distance-to-other-junction filtering is implemented for SAM output.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::out_sj::{Junction, OutSJ};

/// Minimum set of `outSJfilter*` parameters. Extracted from `Parameters`.
#[derive(Debug, Clone)]
pub struct SjFilter {
    /// `--outSJfilterCountUniqueMin`, 4 values indexed by (motif+1)/2.
    pub count_unique_min: [u32; 4],
    /// `--outSJfilterCountTotalMin`.
    pub count_total_min: [u32; 4],
    /// `--outSJfilterOverhangMin`.
    pub overhang_min: [u16; 4],
    /// `--outSJfilterDistToOtherSJmin`.
    pub dist_to_other_sj_min: [u64; 4],
    /// `--outSJfilterIntronMaxVsReadN`.
    pub intron_max_vs_read_n: Vec<u64>,
    /// If `true`, skip outSJfilter chain (corresponds to `outFilterBySJoutStage != 2`).
    pub stage2: bool,
}

impl Default for SjFilter {
    fn default() -> Self {
        Self {
            count_unique_min: [3, 1, 1, 1],
            count_total_min: [3, 1, 1, 1],
            overhang_min: [30, 12, 12, 12],
            dist_to_other_sj_min: [10, 0, 5, 10],
            // C++ parametersDefault: `outSJfilterIntronMaxVsReadN 50000 100000
            // 200000`. A junction supported by 1 read is kept only if its
            // intron is <=50000, with 2 reads <=100000, with 3 <=200000;
            // >=4-read junctions fall back to `alignIntronMax`. Leaving this
            // empty effectively disabled the intron-length cap.
            intron_max_vs_read_n: vec![50_000, 100_000, 200_000],
            stage2: false,
        }
    }
}

/// Merge a set of chunk-level `OutSJ` collectors into `all` (collapsed),
/// then apply filtering and write `SJ.out.tab`.
///
/// `chr_name`/`chr_start` give the chromosome name & start offset for each
/// chr id (parallel to `Genome::chr_name` / `Genome::chr_start`).
pub fn output_sj(
    chunks: &mut [OutSJ],
    filter: &SjFilter,
    chr_name: &[String],
    chr_start: &[u64],
    sj_out_path: impl AsRef<Path>,
) -> std::io::Result<()> {
    let mut all = OutSJ::new();
    for ch in chunks.iter_mut() {
        ch.collapse_sj();
        all.append(ch);
    }
    all.collapse_sj();

    let sj_filter_flags: Vec<bool> = if filter.stage2 {
        vec![true; all.junctions.len()]
    } else {
        compute_filter_flags(&all.junctions, filter)
    };

    let file = File::create(sj_out_path)?;
    let mut writer = BufWriter::new(file);

    for (ii, j) in all.junctions.iter().enumerate() {
        if !sj_filter_flags[ii] && j.annot == 0 {
            continue;
        }
        let chr_idx = chr_idx_for_start(chr_start, j.start);
        OutSJ::write_junction_line(
            &mut writer,
            &chr_name[chr_idx],
            chr_start[chr_idx],
            j,
        )?;
    }
    writer.flush()?;
    Ok(())
}

/// Motif sentinel used by C++ for annotated junctions — beyond any real
/// motif index (motifs are 0..=6), so the `(motif+1)/2` bin computation
/// maps to an out-of-range entry which C++ treats as "no filter". We
/// achieve the same effect by short-circuiting on this sentinel.
const SJ_MOTIF_ANNOTATED: u8 = 7;
/// Port of `outputSJ.cpp` filter chain. Mirrors the C++ two-stage logic:
///
/// 1. **Count / overhang / intron-length filter (`stage1`)** — a junction
///    is considered a candidate only if it passes this test. The C++ code
///    physically drops rejected junctions from `allSJ` before running the
///    distance check; we keep them in `junctions` but track candidacy via
///    a bitmap so the distance pass ignores them.
/// 2. **Donor-neighbor distance filter** — scan junctions sorted by donor
///    position (already the order of `junctions`). For each candidate,
///    min-distance to the nearest *candidate* donor on either side must
///    be >= `distToOtherSJmin[bin]`.
/// 3. **Acceptor-neighbor distance filter** — the same test but with
///    junctions sorted by acceptor position (`start + gap`); both passes
///    AND together.
///
/// Annotated junctions bypass both distance stages (their motif slot is
/// replaced with the sentinel `SJ_MOTIF_ANNOTATED`).
fn compute_filter_flags(junctions: &[Junction], filter: &SjFilter) -> Vec<bool> {
    let n = junctions.len();
    let mut flags = vec![false; n];
    let mut candidate = vec![false; n];

    // Stage 1 — count/overhang/intron-length. Annotated junctions are
    // kept unconditionally.
    for (ii, j) in junctions.iter().enumerate() {
        if j.annot > 0 {
            flags[ii] = true;
            candidate[ii] = true;
            continue;
        }
        let bin = ((j.motif as usize + 1) / 2).min(3);
        let base_pass = (j.count_unique >= filter.count_unique_min[bin]
            || (j.count_multiple + j.count_unique) >= filter.count_total_min[bin])
            && j.overhang_left >= filter.overhang_min[bin]
            && j.overhang_right >= filter.overhang_min[bin]
            && ((j.count_multiple + j.count_unique) as usize
                > filter.intron_max_vs_read_n.len()
                || j.gap as u64
                    <= filter.intron_max_vs_read_n
                        [(j.count_multiple + j.count_unique) as usize - 1]);
        if base_pass {
            candidate[ii] = true;
        }
    }

    // Stage 2 — donor-side distance (junctions are already sorted by
    // (start, gap, strand) from `collapse_sj`). Non-candidate neighbors
    // are invisible to the distance check, matching C++'s `allSJ`.
    // Precompute neighbor-candidate positions by walking the sorted list.
    let mut prev_cand: Vec<Option<u64>> = vec![None; n];
    let mut next_cand: Vec<Option<u64>> = vec![None; n];
    let mut last: Option<u64> = None;
    for ii in 0..n {
        prev_cand[ii] = last;
        if candidate[ii] {
            last = Some(junctions[ii].start);
        }
    }
    let mut next: Option<u64> = None;
    for ii in (0..n).rev() {
        next_cand[ii] = next;
        if candidate[ii] {
            next = Some(junctions[ii].start);
        }
    }

    let donor_pass: Vec<bool> = (0..n)
        .map(|ii| {
            if !candidate[ii] {
                return false;
            }
            if junctions[ii].annot > 0 {
                return true;
            }
            let bin = ((junctions[ii].motif as usize + 1) / 2).min(3);
            let here = junctions[ii].start;
            let x1 = prev_cand[ii].unwrap_or(0);
            let x2 = next_cand[ii].unwrap_or(u64::MAX);
            let min_dist = std::cmp::min(here - x1, x2 - here);
            min_dist >= filter.dist_to_other_sj_min[bin]
        })
        .collect();

    // Stage 3 — acceptor-side distance. Build (acceptor_pos, index,
    // motif-or-annotated-sentinel) tuples for CANDIDATES only, sort by
    // acceptor, then apply the same minDist >= threshold test. Annotated
    // junctions carry the sentinel motif so the `(m+1)/2` bin lookup
    // hits the 4th slot (which holds the non-canonical threshold but is
    // short-circuited below).
    let mut acceptor_sorted: Vec<(u64, usize, u8)> = Vec::new();
    for ii in 0..n {
        if !candidate[ii] {
            continue;
        }
        let motif_slot = if junctions[ii].annot > 0 {
            SJ_MOTIF_ANNOTATED
        } else {
            junctions[ii].motif
        };
        let acceptor = junctions[ii].start + junctions[ii].gap as u64;
        acceptor_sorted.push((acceptor, ii, motif_slot));
    }
    acceptor_sorted.sort_by_key(|&(acc, _, _)| acc);

    // Second-pass flag, indexed by original junction index.
    let mut acceptor_pass: Vec<bool> = vec![false; n];
    for jj in 0..acceptor_sorted.len() {
        let (acc, ii, motif_slot) = acceptor_sorted[jj];
        if motif_slot == SJ_MOTIF_ANNOTATED {
            acceptor_pass[ii] = true;
            continue;
        }
        let bin = ((motif_slot as usize + 1) / 2).min(3);
        let x1 = if jj > 0 { acceptor_sorted[jj - 1].0 } else { 0 };
        let x2 = if jj + 1 < acceptor_sorted.len() {
            acceptor_sorted[jj + 1].0
        } else {
            u64::MAX
        };
        let min_dist = std::cmp::min(acc - x1, x2 - acc);
        acceptor_pass[ii] = min_dist >= filter.dist_to_other_sj_min[bin];
    }

    for ii in 0..n {
        if junctions[ii].annot > 0 {
            flags[ii] = true;
        } else {
            flags[ii] = candidate[ii] && donor_pass[ii] && acceptor_pass[ii];
        }
    }
    flags
}

fn chr_idx_for_start(chr_start: &[u64], start: u64) -> usize {
    match chr_start.binary_search(&start) {
        Ok(idx) => idx,
        Err(idx) => idx.saturating_sub(1),
    }
}
