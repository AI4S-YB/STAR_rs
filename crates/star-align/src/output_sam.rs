//! Simplified port of `ReadAlign::outputTranscriptSAM` +
//! `ReadAlign::outputTranscriptCIGARp` + `ReadAlign::outputAlignments`.
//!
//! Scope for M3 basic SE:
//! - Single-end mapped / unmapped reads only (PE mate-split still TODO).
//! - Standard attribute set: NH, HI, AS, nM, jM, jI, NM, MD.
//! - Writes plain SAM text; BAM is deferred to M7 via noodles.

use std::fmt::Write as _;
use std::io::Write;

use star_core::seq::rev_complement_nucleotides;
use star_core::types::{
    ATTR_AS, ATTR_HI, ATTR_JI, ATTR_JM, ATTR_MD, ATTR_NH, ATTR_NM, ATTR_NM_LOWER, ATTR_XS, EX_G,
    EX_I_FRAG, EX_L, EX_R, SJ_SAM_ANNOTATED_MOTIF_SHIFT,
};
use star_genome::Genome;
use star_params::parameters::Parameters;
use star_sjdb::out_sj::OutSJ;

use crate::read_align::ReadAlign;
use crate::transcript::Transcript;

/// Numeric-to-ASCII base table (matches `Parameters::genomeNumToNT`).
const GENOME_NUM_TO_NT: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

impl ReadAlign {
    /// Port of `ReadAlign::outputTranscriptSAM`.
    ///
    /// Handles both the mapped (1 or 2 mates) and unmapped-in-SAM paths.
    /// The mate split is detected via `canon_sj[iEx] == -3` (mate separator
    /// set by `stitch_pieces`). For PE, each mate is written as a separate
    /// SAM line with the proper FLAG / RNEXT / PNEXT / TLEN fields.
    #[allow(clippy::too_many_arguments)]
    pub fn output_transcript_sam(
        &mut self,
        tr_out: &Transcript,
        n_tr_out: u64,
        i_tr_out: u64,
        mate_chr: u64,
        mate_start: u64,
        mate_strand: i32,
        unmap_type: i32,
        mate_map: &[bool; 2],
        p: &Parameters,
        map_gen: &Genome,
        out: &mut impl Write,
        sam_attr_order: &[u16],
    ) -> std::io::Result<usize> {
        if p.out_sam_mode == "None" {
            return Ok(0);
        }
        let mut bytes = 0usize;
        let read_name_stripped = self
            .read_name
            .strip_prefix('@')
            .or_else(|| self.read_name.strip_prefix('>'))
            .unwrap_or(&self.read_name)
            .to_string();
        let flag_paired = p.read_nmates == 2;
        let n_chr_real = map_gen.n_chr_real;

        // Unmapped path: emit one SAM line per unmapped mate.
        if unmap_type >= 0 {
            for imate in 0..p.read_nmates as usize {
                if mate_map.get(imate).copied().unwrap_or(false) {
                    continue;
                }
                let mut sam_flag: u16 = 0x4;
                if flag_paired {
                    sam_flag |= 0x1;
                    sam_flag |= if imate == 0 { 0x40 } else { 0x80 };
                    let other = 1 - imate;
                    let other_mapped = mate_map.get(other).copied().unwrap_or(false);
                    if other_mapped {
                        if tr_out.str_ as usize != other {
                            sam_flag |= 0x20;
                        }
                    } else {
                        sam_flag |= 0x8;
                    }
                }
                if self.read_filter == b'Y' {
                    sam_flag |= 0x200;
                }
                let other_mapped = if flag_paired {
                    mate_map.get(1 - imate).copied().unwrap_or(false)
                } else {
                    false
                };
                let (rname, rpos) = if other_mapped {
                    (
                        map_gen.chr_name[tr_out.chr as usize].clone(),
                        tr_out.exons[0][EX_G] + 1 - map_gen.chr_start[tr_out.chr as usize],
                    )
                } else {
                    ("*".to_string(), 0u64)
                };
                let qual = if self.read_file_type == 2 {
                    std::str::from_utf8(&self.qual0[imate])
                        .unwrap_or("*")
                        .to_string()
                } else {
                    "*".to_string()
                };
                let seq = std::str::from_utf8(&self.read0[imate])
                    .unwrap_or("*")
                    .to_string();
                let line = format!(
                    "{}\t{}\t*\t0\t0\t*\t{}\t{}\t0\t{}\t{}\tNH:i:0\tHI:i:0\tAS:i:{}\tnM:i:{}\tuT:A:{}\n",
                    read_name_stripped,
                    sam_flag,
                    rname,
                    rpos,
                    seq,
                    qual,
                    tr_out.max_score,
                    tr_out.n_mm,
                    unmap_type,
                );
                out.write_all(line.as_bytes())?;
                bytes += line.len();
            }
            return Ok(bytes);
        }

        // Mapped path: split mates by the -3 separator.
        let mut i_ex_mate: u64 = tr_out.n_exons.saturating_sub(1);
        let mut n_mates: u64 = 1;
        for i in 0..tr_out.n_exons.saturating_sub(1) {
            if tr_out.canon_sj[i as usize] == -3 {
                i_ex_mate = i;
                n_mates = 2;
                break;
            }
        }

        // samFlagCommon
        let mut sam_flag_common: u16 = 0;
        if flag_paired {
            sam_flag_common = 0x1;
            if i_ex_mate == tr_out.n_exons - 1 {
                if mate_chr > n_chr_real {
                    sam_flag_common |= 0x8;
                }
            } else {
                // properly paired check: alignEndsProtrude.concordantPair OR coord check.
                let a0g = tr_out.exons[0][EX_G];
                let a0r = tr_out.exons[0][EX_R];
                let m1g = tr_out.exons[i_ex_mate as usize + 1][EX_G];
                let a_mid_g = tr_out.exons[i_ex_mate as usize][EX_G];
                let a_mid_l = tr_out.exons[i_ex_mate as usize][EX_L];
                let last_g = tr_out.exons[tr_out.n_exons as usize - 1][EX_G];
                let last_r = tr_out.exons[tr_out.n_exons as usize - 1][EX_R];
                let concordant = p.align_ends_protrude.concordant_pair
                    || (a0g <= m1g + a0r && a_mid_g + a_mid_l <= last_g + self.l_read - last_r);
                if concordant {
                    sam_flag_common |= 0x2;
                }
            }
        }
        if self.read_filter == b'Y' {
            sam_flag_common |= 0x200;
        }

        let str_ = tr_out.str_;
        let left_mate = if flag_paired { str_ } else { 0 };

        let mapq = if n_tr_out >= 5 {
            0
        } else if n_tr_out >= 3 {
            1
        } else if n_tr_out == 2 {
            3
        } else {
            p.out_sam_mapq_unique as i32
        };

        let (tag_nm, tag_md) =
            if sam_attr_order.contains(&ATTR_NM) || sam_attr_order.contains(&ATTR_MD) {
                // NM/MD are per whole transcript in C++; we compute once per mate below.
                (0u64, String::new())
            } else {
                (0u64, String::new())
            };
        let _ = (tag_nm, tag_md);

        for imate in 0..n_mates {
            let i_ex1 = if imate == 0 { 0 } else { i_ex_mate + 1 };
            let i_ex2 = if imate == 0 {
                i_ex_mate
            } else {
                tr_out.n_exons - 1
            };
            let mate = tr_out.exons[i_ex1 as usize][EX_I_FRAG];

            let mut sam_flag = sam_flag_common;
            if mate == 0 {
                sam_flag |= (str_ as u16) * 0x10;
                if n_mates == 2 {
                    sam_flag |= ((1 - str_) as u16) * 0x20;
                }
            } else {
                sam_flag |= ((1 - str_) as u16) * 0x10;
                if n_mates == 2 {
                    sam_flag |= (str_ as u16) * 0x20;
                }
            }
            if flag_paired {
                sam_flag |= if mate == 0 { 0x40 } else { 0x80 };
                if n_mates == 1 && mate_strand == 1 {
                    sam_flag |= 0x20;
                }
            }
            if !tr_out.primary_flag {
                sam_flag |= 0x100;
            }

            let (cigar, sj_motif_str, sj_intron_str) =
                build_cigar_and_sj_range(self, tr_out, i_ex1, i_ex2, mate, left_mate, map_gen);

            let (seq_str, qual_str) = mate_seq_and_qual_imate(self, str_, mate);

            let flag_final = (sam_flag & p.out_sam_flag_and as u16) | p.out_sam_flag_or as u16;
            let chr_name = &map_gen.chr_name[tr_out.chr as usize];
            let pos_1based =
                tr_out.exons[i_ex1 as usize][EX_G] + 1 - map_gen.chr_start[tr_out.chr as usize];

            write!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}",
                read_name_stripped, flag_final, chr_name, pos_1based, mapq, cigar
            )?;

            if n_mates > 1 {
                // RNEXT == "=", PNEXT == mate pos, TLEN == fragment length.
                let other_ex = if imate == 0 { i_ex_mate + 1 } else { 0 };
                let pnext = tr_out.exons[other_ex as usize][EX_G] + 1
                    - map_gen.chr_start[tr_out.chr as usize];
                let tlen_val = tr_out.exons[tr_out.n_exons as usize - 1][EX_G]
                    + tr_out.exons[tr_out.n_exons as usize - 1][EX_L]
                    - tr_out.exons[0][EX_G];
                let sign = if imate == 0 { "" } else { "-" };
                write!(out, "\t=\t{}\t{}{}", pnext, sign, tlen_val)?;
            } else if mate_chr < n_chr_real {
                let rnext = &map_gen.chr_name[mate_chr as usize];
                let pnext = mate_start + 1 - map_gen.chr_start[mate_chr as usize];
                write!(out, "\t{}\t{}\t0", rnext, pnext)?;
            } else {
                write!(out, "\t*\t0\t0")?;
            }

            write!(out, "\t{}", seq_str)?;
            if self.read_file_type == 2 && p.out_sam_mode != "NoQS" {
                write!(out, "\t{}", qual_str)?;
            } else {
                write!(out, "\t*")?;
            }

            let (tag_nm, tag_md) =
                if sam_attr_order.contains(&ATTR_NM) || sam_attr_order.contains(&ATTR_MD) {
                    compute_nm_md(self, tr_out, i_ex1, i_ex2, map_gen)
                } else {
                    (0u64, String::new())
                };

            let mut buf = String::new();
            for attr in sam_attr_order {
                match *attr {
                    ATTR_NH => write!(buf, "\tNH:i:{}", n_tr_out).unwrap(),
                    ATTR_HI => write!(buf, "\tHI:i:{}", i_tr_out + 1).unwrap(),
                    ATTR_AS => write!(buf, "\tAS:i:{}", tr_out.max_score).unwrap(),
                    ATTR_NM_LOWER => write!(buf, "\tnM:i:{}", tr_out.n_mm).unwrap(),
                    ATTR_JM => write!(buf, "\tjM:B:c{}", sj_motif_str).unwrap(),
                    ATTR_JI => write!(buf, "\tjI:B:i{}", sj_intron_str).unwrap(),
                    ATTR_XS => match tr_out.sj_motif_strand {
                        1 => buf.push_str("\tXS:A:+"),
                        2 => buf.push_str("\tXS:A:-"),
                        _ => {}
                    },
                    ATTR_NM => write!(buf, "\tNM:i:{}", tag_nm).unwrap(),
                    ATTR_MD => write!(buf, "\tMD:Z:{}", tag_md).unwrap(),
                    _ => {}
                }
            }
            out.write_all(buf.as_bytes())?;
            out.write_all(b"\n")?;
            bytes += buf.len() + 1;
        }
        Ok(bytes)
    }

    /// 1:1 port of `ReadAlign::outputTranscriptCIGARp`
    /// (ReadAlign_outputTranscriptCIGARp.cpp) — produces the "CIGAR + p"
    /// string used for `Chimeric.out.junction` columns 12/14.
    pub fn output_transcript_cigar_p(&self, tr_out: &Transcript, p: &Parameters) -> String {
        let left_mate = if p.read_files_in.len() > 1 {
            tr_out.str_ as usize
        } else {
            0
        };
        let lr_orig_left = self.read_length_original[left_mate];
        let lr_orig_pair = self.read_length_pair_original;
        let mut cigar = String::new();
        let e0_r = tr_out.exons[0][EX_R];
        let sub = if e0_r < lr_orig_left {
            0
        } else {
            lr_orig_left + 1
        };
        let trim_l = e0_r - sub;
        if trim_l > 0 {
            write!(cigar, "{}S", trim_l).unwrap();
        }
        for ii in 0..tr_out.n_exons as usize {
            if ii > 0 {
                let prev = &tr_out.exons[ii - 1];
                let cur = &tr_out.exons[ii];
                let gap_g = cur[EX_G] as i64 - (prev[EX_G] + prev[EX_L]) as i64;
                if gap_g >= 0 {
                    if tr_out.canon_sj[ii - 1] == -3 {
                        let s1 = lr_orig_left - (prev[EX_R] + prev[EX_L]);
                        let s2 = cur[EX_R] - (lr_orig_left + 1);
                        if s1 > 0 {
                            write!(cigar, "{}S", s1).unwrap();
                        }
                        write!(cigar, "{}p", gap_g).unwrap();
                        if s2 > 0 {
                            write!(cigar, "{}S", s2).unwrap();
                        }
                    } else {
                        let gap_r = cur[EX_R] - prev[EX_R] - prev[EX_L];
                        if gap_r > 0 {
                            write!(cigar, "{}I", gap_r).unwrap();
                        }
                        if tr_out.canon_sj[ii - 1] >= 0 || tr_out.sj_annot[ii - 1] == 1 {
                            write!(cigar, "{}N", gap_g).unwrap();
                        } else if gap_g > 0 {
                            write!(cigar, "{}D", gap_g).unwrap();
                        }
                    }
                } else {
                    let ov = (prev[EX_G] + prev[EX_L]) as i64 - cur[EX_G] as i64;
                    write!(cigar, "-{}p", ov).unwrap();
                }
            }
            write!(cigar, "{}M", tr_out.exons[ii][EX_L]).unwrap();
        }
        let last = tr_out.n_exons as usize - 1;
        let lr_end_base = if tr_out.exons[last][EX_R] < lr_orig_left {
            lr_orig_left
        } else {
            lr_orig_pair
        };
        let trim_r = lr_end_base - tr_out.exons[last][EX_R] - tr_out.exons[last][EX_L];
        if trim_r > 0 {
            write!(cigar, "{}S", trim_r).unwrap();
        }
        cigar
    }

    /// Port of `ReadAlign::outputAlignments` — basic SE path; writes SAM +
    /// records SJs.
    pub fn output_alignments(
        &mut self,
        p: &Parameters,
        map_gen: &Genome,
        sam_out: &mut impl Write,
        chunk_out_sj: &mut OutSJ,
        sam_attr_order: &[u16],
    ) -> anyhow::Result<()> {
        // Record SJs for passing aligns.
        if self.unmap_type < 0 && self.n_tr > 0 {
            let sj_start_n = chunk_out_sj.n();
            for tr in &self.tr_mult_array {
                let exons: Vec<(u64, u64)> = tr
                    .exons
                    .iter()
                    .take(tr.n_exons as usize)
                    .map(|e| (e[EX_G], e[EX_L]))
                    .collect();
                chunk_out_sj.record_transcript_sjs(
                    &exons,
                    &tr.canon_sj[..tr.n_exons as usize],
                    &tr.sj_annot[..tr.n_exons as usize],
                    self.n_tr,
                    sj_start_n,
                );
            }
        }

        // Write SAM.
        if self.unmap_type < 0 && self.n_tr > 0 {
            let tr_mult = std::mem::take(&mut self.tr_mult_array);
            for (i_tr, tr) in tr_mult.iter().enumerate() {
                let mate_map = [true, false];
                self.output_transcript_sam(
                    tr,
                    self.n_tr,
                    i_tr as u64,
                    u64::MAX,
                    0,
                    0,
                    -1,
                    &mate_map,
                    p,
                    map_gen,
                    sam_out,
                    sam_attr_order,
                )?;
            }
            self.tr_mult_array = tr_mult;
        } else if self.unmap_type >= 0 && p.out_sam_unmapped == "Within" {
            // Write a dummy SAM line for unmapped read.
            let dummy = crate::transcript::Transcript::new();
            let mate_map = [false, false];
            self.output_transcript_sam(
                &dummy,
                0,
                0,
                u64::MAX,
                0,
                0,
                self.unmap_type,
                &mate_map,
                p,
                map_gen,
                sam_out,
                sam_attr_order,
            )?;
        }

        Ok(())
    }
}

/// Mate-aware port of the CIGAR-building loop from
/// `ReadAlign_outputTranscriptSAM.cpp`. The caller supplies the full `[i_ex1,
/// i_ex2]` range for one mate plus `left_mate` (the mate that is on the left
/// of the fragment, i.e. `str_` when paired).
fn build_cigar_and_sj_range(
    ra: &ReadAlign,
    tr_out: &Transcript,
    i_ex1: u64,
    i_ex2: u64,
    mate: u64,
    left_mate: u64,
    map_gen: &Genome,
) -> (String, String, String) {
    let str_ = tr_out.str_;
    let trim_l = if str_ == 0 && mate == 0 {
        ra.clip_mates[mate as usize][0].clipped_n
    } else if str_ == 0 && mate == 1 {
        ra.clip_mates[mate as usize][1].clipped_n
    } else if str_ == 1 && mate == 0 {
        ra.clip_mates[mate as usize][1].clipped_n
    } else {
        ra.clip_mates[mate as usize][0].clipped_n
    };
    let lr_left = ra.read_length[left_mate as usize];
    let ex1_r = tr_out.exons[i_ex1 as usize][EX_R];
    let sub = if ex1_r < lr_left { 0 } else { lr_left + 1 };
    let trim_l1 = trim_l + ex1_r - sub;

    let mut cigar = String::new();
    let mut motif = String::new();
    let mut intron = String::new();
    if trim_l1 > 0 {
        write!(cigar, "{}S", trim_l1).unwrap();
    }
    for ii in i_ex1..=i_ex2 {
        let iiu = ii as usize;
        if ii > i_ex1 {
            let prev = &tr_out.exons[iiu - 1];
            let cur = &tr_out.exons[iiu];
            let gap_g = cur[EX_G] - (prev[EX_G] + prev[EX_L]);
            let gap_r = cur[EX_R] - prev[EX_R] - prev[EX_L];
            if gap_r > 0 {
                write!(cigar, "{}I", gap_r).unwrap();
            }
            if tr_out.canon_sj[iiu - 1] >= 0 || tr_out.sj_annot[iiu - 1] == 1 {
                write!(cigar, "{}N", gap_g).unwrap();
                let cj = tr_out.canon_sj[iiu - 1];
                let motif_val = if tr_out.sj_annot[iiu - 1] == 0 {
                    cj
                } else {
                    cj + SJ_SAM_ANNOTATED_MOTIF_SHIFT as i32
                };
                write!(motif, ",{}", motif_val).unwrap();
                let istart = prev[EX_G] + prev[EX_L] + 1 - map_gen.chr_start[tr_out.chr as usize];
                let iend = cur[EX_G] - map_gen.chr_start[tr_out.chr as usize];
                write!(intron, ",{},{}", istart, iend).unwrap();
            } else if gap_g > 0 {
                write!(cigar, "{}D", gap_g).unwrap();
            }
        }
        write!(cigar, "{}M", tr_out.exons[iiu][EX_L]).unwrap();
    }

    let lr_orig_left = ra.read_length_original[left_mate as usize];
    let lr_orig_mate = ra.read_length_original[mate as usize];
    let ex2_r = tr_out.exons[i_ex2 as usize][EX_R];
    let ex2_l = tr_out.exons[i_ex2 as usize][EX_L];
    let base = if ex1_r < lr_left {
        lr_orig_left
    } else {
        lr_left + 1 + lr_orig_mate
    };
    let trim_r1 = base.saturating_sub(ex2_r + ex2_l + trim_l);
    if trim_r1 > 0 {
        write!(cigar, "{}S", trim_r1).unwrap();
    }

    if motif.is_empty() {
        motif.push_str(",-1");
        intron.push_str(",-1");
    }

    (cigar, motif, intron)
}

/// Mate-aware port of the seq/qual output from
/// `ReadAlign_outputTranscriptSAM.cpp`. If `Mate == Str` the read is in the
/// correct orientation; otherwise reverse-complement it.
fn mate_seq_and_qual_imate(ra: &ReadAlign, str_: u64, mate: u64) -> (String, String) {
    let seq0 = &ra.read0[mate as usize];
    let qual0 = &ra.qual0[mate as usize];
    if mate == str_ {
        (
            String::from_utf8_lossy(seq0).into_owned(),
            String::from_utf8_lossy(qual0).into_owned(),
        )
    } else {
        let mut seq_rc = vec![0u8; seq0.len()];
        rev_complement_nucleotides(seq0, &mut seq_rc);
        let mut qual_rev = qual0.clone();
        qual_rev.reverse();
        (
            String::from_utf8_lossy(&seq_rc).into_owned(),
            String::from_utf8_lossy(&qual_rev).into_owned(),
        )
    }
}

fn compute_nm_md(
    ra: &ReadAlign,
    tr_out: &Transcript,
    i_ex1: u64,
    i_ex2: u64,
    map_gen: &Genome,
) -> (u64, String) {
    let r = if tr_out.ro_str == 0 {
        &ra.read1[0]
    } else {
        &ra.read1[2]
    };
    let mut tag_nm: u64 = 0;
    let mut tag_md = String::new();
    let mut match_n: u64 = 0;
    for iex in i_ex1..=i_ex2 {
        let e = &tr_out.exons[iex as usize];
        for ii in 0..e[EX_L] {
            let r1 = r[(ii + e[EX_R]) as usize];
            let g1 = map_gen.g[(ii + e[EX_G]) as usize + star_genome::load::LOAD_L];
            if r1 != g1 || r1 == 4 || g1 == 4 {
                tag_nm += 1;
                write!(tag_md, "{}", match_n).unwrap();
                tag_md.push(GENOME_NUM_TO_NT[g1.min(4) as usize] as char);
                match_n = 0;
            } else {
                match_n += 1;
            }
        }
        if iex < i_ex2 {
            let cj = tr_out.canon_sj[iex as usize];
            if cj == -1 {
                let next = &tr_out.exons[iex as usize + 1];
                tag_nm += next[EX_G] - (e[EX_G] + e[EX_L]);
                write!(tag_md, "{}^", match_n).unwrap();
                for ii in e[EX_G] + e[EX_L]..next[EX_G] {
                    tag_md.push(
                        GENOME_NUM_TO_NT
                            [map_gen.g[ii as usize + star_genome::load::LOAD_L].min(4) as usize]
                            as char,
                    );
                }
                match_n = 0;
            } else if cj == -2 {
                let next = &tr_out.exons[iex as usize + 1];
                tag_nm += next[EX_R] - e[EX_R] - e[EX_L];
            }
        }
    }
    write!(tag_md, "{}", match_n).unwrap();
    (tag_nm, tag_md)
}
