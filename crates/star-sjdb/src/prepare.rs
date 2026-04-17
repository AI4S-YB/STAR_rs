//! 1:1 port of `sjdbPrepare.cpp` (STAR/source/sjdbPrepare.cpp:1-225).
//!
//! Input: `SjdbLoci` (accumulated junctions), `Genome` (for G/chr lookup),
//! a `P` for logging/chr-bin params, an output directory for
//! `sjdbInfo.txt` and `sjdbList.out.tab`, and an output `Gsj` byte buffer
//! for the concatenated sjdb sequences (size `2*sjdbLength*nSJloci + 1`).
//!
//! The port preserves the two-stage de-duplication ordering:
//! 1. Sort by `(sjdbS+shift1, sjdbE+shift1, ii)` — mirrors `qsort` with
//!    `funCompareUint2`. `shift1` separates strands into `[0, nG, 2nG)`
//!    buckets (`+`, `-`, `.`).
//! 2. Walk in sorted order, keep the "winner" per intron coordinate based
//!    on priority and canonical/left-most tie-break rules.
//! 3. Second sort on the winners using their *unshifted* canonical loci.
//! 4. Merge opposite-strand duplicates with priority/strand/motif tie-break
//!    rules.
//!
//! After de-dup, the function:
//! - Fills `mapGen.sjdbStart/End/Motif/ShiftLeft/ShiftRight/Strand` +
//!   `sjDstart/sjAstart`.
//! - Writes `outDir/sjdbInfo.txt` (header + one line per junction).
//! - Writes `outDir/sjdbList.out.tab` (chr, start, end, strand char).
//! - Fills `Gsj` with the concatenated donor+acceptor sequences of each
//!   junction (each padded with a trailing `GENOME_spacingChar`).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use star_core::types::GENOME_SPACING_CHAR;
use star_genome::genome::Genome;
use star_genome::load::LOAD_L;
use star_params::parameters::Parameters;

use crate::sjdb_class::SjdbLoci;

/// Port of `sjdbPrepare(SjdbClass&, Parameters&, uint, string, Genome&,
/// char*)`.
///
/// - `n_genome_real`: `mapGen.chrStart[mapGen.nChrReal]` (i.e. the last
///   chromosome-boundary value — used only by the repeat-shift search).
/// - `out_dir`: path to write `sjdbInfo.txt` + `sjdbList.out.tab`.
/// - On return the caller gets a `Vec<u8>` with the concatenated sjdb
///   donor+acceptor sequences.
pub fn sjdb_prepare(
    sjdb_loci: &SjdbLoci,
    p: &Parameters,
    n_genome_real: u64,
    out_dir: impl AsRef<Path>,
    map_gen: &mut Genome,
) -> anyhow::Result<Vec<u8>> {
    let n_in = sjdb_loci.len();
    let mut sjdb_s = vec![0u64; n_in];
    let mut sjdb_e = vec![0u64; n_in];
    let mut sjdb_motif = vec![0u8; n_in];
    let mut sjdb_shift_left = vec![0u8; n_in];
    let mut sjdb_shift_right = vec![0u8; n_in];

    // `G` is the numeric genome view (past the LOAD_L pre-padding added
    // by `genome_load`).
    let g_start = LOAD_L;
    // G access lambda
    let g_at = |idx: u64| -> u8 { map_gen.g[g_start + idx as usize] };

    // ------------------------------------------------------------------
    // Step 1: resolve chromosomes, compute motif + shift_left/shift_right
    // ------------------------------------------------------------------
    let mut chr_old: &str = "";
    let mut i_chr: usize = 0;
    for ii in 0..n_in {
        if chr_old != sjdb_loci.chr[ii] {
            i_chr = map_gen.n_chr_real as usize;
            for c in 0..map_gen.n_chr_real as usize {
                if sjdb_loci.chr[ii] == map_gen.chr_name[c] {
                    i_chr = c;
                    break;
                }
            }
            if i_chr >= map_gen.n_chr_real as usize {
                anyhow::bail!(
                    "EXITING because of FATAL error, the sjdb chromosome {} is not found among the genomic chromosomes\n\
                     SOLUTION: fix your file(s) --sjdbFileChrStartEnd or --sjdbGTFfile, offending junction: {}\t{}\t{}",
                    sjdb_loci.chr[ii], sjdb_loci.chr[ii], sjdb_loci.start[ii], sjdb_loci.end[ii]
                );
            }
            chr_old = &sjdb_loci.chr[ii];
        }

        sjdb_s[ii] = sjdb_loci.start[ii] + map_gen.chr_start[i_chr] - 1;
        sjdb_e[ii] = sjdb_loci.end[ii] + map_gen.chr_start[i_chr] - 1;

        let s = sjdb_s[ii];
        let e = sjdb_e[ii];
        sjdb_motif[ii] = match (g_at(s), g_at(s + 1), g_at(e - 1), g_at(e)) {
            (2, 3, 0, 2) => 1, // GTAG
            (1, 3, 0, 1) => 2, // CTAC
            (2, 1, 0, 2) => 3, // GCAG
            (1, 3, 2, 1) => 4, // CTGC
            (0, 3, 0, 1) => 5, // ATAC
            (2, 3, 0, 3) => 6, // GTAT
            _ => 0,
        };

        // Repeat length: go back around the junction to find max shift left.
        let mut jj_l: u64 = 0;
        while jj_l <= s - 1
            && g_at(s - 1 - jj_l) == g_at(e - jj_l)
            && g_at(s - 1 - jj_l) < 4
            && jj_l < 255
        {
            jj_l += 1;
        }
        sjdb_shift_left[ii] = jj_l as u8;

        let mut jj_r: u64 = 0;
        while s + jj_r < n_genome_real
            && g_at(s + jj_r) == g_at(e + 1 + jj_r)
            && g_at(s + jj_r) < 4
            && jj_r < 255
        {
            jj_r += 1;
        }
        sjdb_shift_right[ii] = jj_r as u8;

        if jj_r == 255 || jj_l == 255 {
            eprintln!(
                "WARNING: long repeat for junction # {} : {} {} {}; left shift = {}; right shift = {}",
                ii + 1,
                sjdb_loci.chr[ii],
                s - map_gen.chr_start[i_chr] + 1,
                e - map_gen.chr_start[i_chr] + 1,
                jj_l,
                jj_r
            );
        }

        sjdb_s[ii] -= sjdb_shift_left[ii] as u64;
        sjdb_e[ii] -= sjdb_shift_left[ii] as u64;
    }

    // ------------------------------------------------------------------
    // Step 2: first sort (stranded) + priority/canonical tie-break dedupe
    // ------------------------------------------------------------------
    let n_g = n_genome_real;
    let mut sort_a: Vec<(u64, u64, u64)> = (0..n_in as u64)
        .map(|ii| {
            let shift1 = match sjdb_loci.str_[ii as usize] {
                b'+' => 0u64,
                b'-' => n_g,
                _ => 2 * n_g,
            };
            (
                sjdb_s[ii as usize] + shift1,
                sjdb_e[ii as usize] + shift1,
                ii,
            )
        })
        .collect();
    sort_a.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    let mut winners: Vec<u64> = Vec::with_capacity(n_in);
    for &(_, _, isj) in &sort_a {
        let isj = isj as usize;
        if winners.is_empty() {
            winners.push(isj as u64);
            continue;
        }
        let last_idx = *winners.last().unwrap() as usize;
        if sjdb_s[isj] != sjdb_s[last_idx] || sjdb_e[isj] != sjdb_e[last_idx] {
            winners.push(isj as u64);
            continue;
        }
        // Same intron coordinates — apply priority & canonical tie-break.
        let pri_new = sjdb_loci.priority[isj];
        let pri_old = sjdb_loci.priority[last_idx];
        if pri_new < pri_old {
            // keep old
        } else if pri_new > pri_old {
            *winners.last_mut().unwrap() = isj as u64;
        } else {
            let new_is_canon = sjdb_motif[isj] > 0;
            let old_is_canon = sjdb_motif[last_idx] > 0;
            let new_wins = (new_is_canon && !old_is_canon)
                || (new_is_canon == old_is_canon && sjdb_shift_left[isj] < sjdb_shift_left[last_idx]);
            if new_wins {
                *winners.last_mut().unwrap() = isj as u64;
            }
        }
    }

    // ------------------------------------------------------------------
    // Step 3: second sort using unshifted canonical loci
    // ------------------------------------------------------------------
    let mut sort_b: Vec<(u64, u64, u64)> = winners
        .iter()
        .map(|&isj| {
            let isj_us = isj as usize;
            let shift = if sjdb_motif[isj_us] == 0 {
                0
            } else {
                sjdb_shift_left[isj_us] as u64
            };
            (sjdb_s[isj_us] + shift, sjdb_e[isj_us] + shift, isj)
        })
        .collect();
    sort_b.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    // ------------------------------------------------------------------
    // Step 4: merge opposite-strand duplicates
    // ------------------------------------------------------------------
    let n_win = sort_b.len();
    let mut out_start: Vec<u64> = Vec::with_capacity(n_win);
    let mut out_end: Vec<u64> = Vec::with_capacity(n_win);
    let mut out_motif: Vec<u8> = Vec::with_capacity(n_win);
    let mut out_shift_l: Vec<u8> = Vec::with_capacity(n_win);
    let mut out_shift_r: Vec<u8> = Vec::with_capacity(n_win);
    let mut out_strand: Vec<u8> = Vec::with_capacity(n_win);

    for ii in 0..n_win {
        let isj = sort_b[ii].2 as usize;
        let s_here = sort_b[ii].0;
        let e_here = sort_b[ii].1;

        if !out_start.is_empty() {
            let last = out_start.len() - 1;
            if out_start[last] == s_here && out_end[last] == e_here {
                let isj0 = sort_b[ii - 1].2 as usize;
                let pri_new = sjdb_loci.priority[isj];
                let pri_old = sjdb_loci.priority[isj0];
                if pri_new < pri_old {
                    continue;
                } else if pri_new > pri_old {
                    out_start.truncate(last);
                    out_end.truncate(last);
                    out_motif.truncate(last);
                    out_shift_l.truncate(last);
                    out_shift_r.truncate(last);
                    out_strand.truncate(last);
                } else if out_strand[last] > 0 && sjdb_loci.str_[isj] == b'.' {
                    continue;
                } else if out_strand[last] == 0 && sjdb_loci.str_[isj] != b'.' {
                    out_start.truncate(last);
                    out_end.truncate(last);
                    out_motif.truncate(last);
                    out_shift_l.truncate(last);
                    out_shift_r.truncate(last);
                    out_strand.truncate(last);
                } else if out_motif[last] == 0 && sjdb_motif[isj] == 0 {
                    // both non-canonical -> mark undefined strand, keep old
                    out_strand[last] = 0;
                    continue;
                } else if (out_motif[last] > 0 && sjdb_motif[isj] == 0)
                    || (out_motif[last] % 2 == (2 - out_strand[last]))
                {
                    continue;
                } else {
                    out_start.truncate(last);
                    out_end.truncate(last);
                    out_motif.truncate(last);
                    out_shift_l.truncate(last);
                    out_shift_r.truncate(last);
                    out_strand.truncate(last);
                }
            }
        }

        out_start.push(s_here);
        out_end.push(e_here);
        out_motif.push(sjdb_motif[isj]);
        out_shift_l.push(sjdb_shift_left[isj]);
        out_shift_r.push(sjdb_shift_right[isj]);
        out_strand.push(match sjdb_loci.str_[isj] {
            b'+' => 1,
            b'-' => 2,
            _ => {
                if sjdb_motif[isj] == 0 {
                    0
                } else {
                    2 - sjdb_motif[isj] % 2
                }
            }
        });
    }

    // ------------------------------------------------------------------
    // Step 5: commit to Genome, build Gsj, write sjdbInfo.txt + sjdbList.out.tab
    // ------------------------------------------------------------------
    let nsj = out_start.len() as u64;
    map_gen.sjdb_n = nsj;
    map_gen.sjdb_start = out_start.clone();
    map_gen.sjdb_end = out_end.clone();
    map_gen.sjdb_motif = out_motif.clone();
    map_gen.sjdb_shift_left = out_shift_l.clone();
    map_gen.sjdb_shift_right = out_shift_r.clone();
    map_gen.sjdb_strand = out_strand.clone();
    map_gen.sj_d_start = vec![0u64; nsj as usize];
    map_gen.sj_a_start = vec![0u64; nsj as usize];

    let sjdb_overhang = map_gen.sjdb_overhang;
    let sjdb_length = map_gen.sjdb_length;
    if sjdb_length == 0 {
        // `sjdbLength = 2 * sjdbOverhang + 1` — set at genome generation.
        // For runtime-inserted junctions the caller should set it before
        // calling sjdb_prepare.
        anyhow::bail!("sjdb_prepare called with genome.sjdb_length == 0");
    }
    let gsj_len = (2 * sjdb_length * nsj + 1) as usize;
    let mut gsj = vec![0u8; gsj_len];

    let out_dir = out_dir.as_ref();
    std::fs::create_dir_all(out_dir)?;
    let mut sjdb_info =
        BufWriter::new(File::create(out_dir.join("sjdbInfo.txt")).map_err(|e| {
            anyhow::anyhow!(
                "EXITING because of fatal ERROR: could not open {}/sjdbInfo.txt: {e}",
                out_dir.display()
            )
        })?);
    let mut sjdb_list =
        BufWriter::new(File::create(out_dir.join("sjdbList.out.tab")).map_err(|e| {
            anyhow::anyhow!(
                "EXITING because of fatal ERROR: could not open {}/sjdbList.out.tab: {e}",
                out_dir.display()
            )
        })?);

    writeln!(sjdb_info, "{}\t{}", nsj, sjdb_overhang)?;
    let strand_char = [b'.', b'+', b'-'];

    let mut sj_g_start: u64 = 0;
    for ii in 0..nsj as usize {
        map_gen.sj_d_start[ii] = map_gen.sjdb_start[ii] - sjdb_overhang;
        map_gen.sj_a_start[ii] = map_gen.sjdb_end[ii] + 1;
        if map_gen.sjdb_motif[ii] == 0 {
            map_gen.sj_d_start[ii] += map_gen.sjdb_shift_left[ii] as u64;
            map_gen.sj_a_start[ii] += map_gen.sjdb_shift_left[ii] as u64;
        }
        // Gsj: donor + acceptor regions, each `sjdbOverhang` bytes, then one
        // spacer char.
        let d_src = (LOAD_L as u64 + map_gen.sj_d_start[ii]) as usize;
        let a_src = (LOAD_L as u64 + map_gen.sj_a_start[ii]) as usize;
        let dst_d = sj_g_start as usize;
        let dst_a = (sj_g_start + sjdb_overhang) as usize;
        gsj[dst_d..dst_d + sjdb_overhang as usize]
            .copy_from_slice(&map_gen.g[d_src..d_src + sjdb_overhang as usize]);
        gsj[dst_a..dst_a + sjdb_overhang as usize]
            .copy_from_slice(&map_gen.g[a_src..a_src + sjdb_overhang as usize]);
        sj_g_start += sjdb_length;
        gsj[sj_g_start as usize - 1] = GENOME_SPACING_CHAR;

        writeln!(
            sjdb_info,
            "{}\t{}\t{}\t{}\t{}\t{}",
            map_gen.sjdb_start[ii],
            map_gen.sjdb_end[ii],
            map_gen.sjdb_motif[ii] as i32,
            map_gen.sjdb_shift_left[ii] as i32,
            map_gen.sjdb_shift_right[ii] as i32,
            map_gen.sjdb_strand[ii] as i32
        )?;

        // sjdbList line — convert to 1-based and shift by chrStart.
        let chr1_idx = (map_gen.sjdb_start[ii] >> p.p_ge.g_chr_bin_nbits) as usize;
        let chr1 = map_gen.chr_bin[chr1_idx] as usize;
        let shift = if map_gen.sjdb_motif[ii] > 0 {
            0
        } else {
            map_gen.sjdb_shift_left[ii] as u64
        };
        writeln!(
            sjdb_list,
            "{}\t{}\t{}\t{}",
            map_gen.chr_name[chr1],
            map_gen.sjdb_start[ii] - map_gen.chr_start[chr1] + 1 + shift,
            map_gen.sjdb_end[ii] - map_gen.chr_start[chr1] + 1 + shift,
            strand_char[map_gen.sjdb_strand[ii] as usize] as char
        )?;
    }
    sjdb_info.flush()?;
    sjdb_list.flush()?;

    Ok(gsj)
}
