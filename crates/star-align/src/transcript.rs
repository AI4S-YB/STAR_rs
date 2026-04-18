//! 1:1 port of `class Transcript` (Transcript.h).
//!
//! Field names match the C++ exactly (with `snake_case` conversion). Array
//! sizes are `MAX_N_EXONS` from star-core; the `long-reads` feature flag
//! scales that to 1000 matching `COMPILE_FOR_LONG_READS`.

use star_core::types::{IntScore, MAX_N_EXONS};
use star_core::{UInt, EX_G, EX_L, EX_R, EX_SIZE};

#[derive(Clone, Debug)]
pub struct Transcript {
    /// `exons[MAX_N_EXONS][EX_SIZE]`: r-start, g-start, length, ...
    pub exons: Vec<[UInt; EX_SIZE]>,
    /// `shiftSJ[MAX_N_EXONS][2]`.
    pub shift_sj: Vec<[UInt; 2]>,
    /// `canonSJ[MAX_N_EXONS]`.
    pub canon_sj: Vec<i32>,
    /// `sjAnnot[MAX_N_EXONS]`.
    pub sj_annot: Vec<u8>,
    /// `sjStr[MAX_N_EXONS]`.
    pub sj_str: Vec<u8>,
    pub intron_motifs: [UInt; 3],
    pub sj_motif_strand: u8,
    pub sj_yes: bool,

    pub n_exons: UInt,

    pub l_read: UInt,
    pub read_length_pair_original: UInt,
    pub i_read: UInt,
    pub read_n_mates: UInt,

    pub i_frag: i32,

    pub r_start: UInt,
    pub ro_start: UInt,
    pub r_length: UInt,
    pub g_start: UInt,
    pub g_length: UInt,
    pub c_start: UInt,
    pub chr: UInt,
    pub str_: UInt,
    pub ro_str: UInt,
    pub haplo_type: u32,

    pub primary_flag: bool,

    pub n_match: UInt,
    pub n_mm: UInt,
    pub mapped_length: UInt,
    pub extend_l: UInt,
    pub max_score: IntScore,

    pub n_gap: UInt,
    pub l_gap: UInt,
    pub n_del: UInt,
    pub n_ins: UInt,
    pub l_del: UInt,
    pub l_ins: UInt,

    pub n_unique: UInt,
    pub n_anchor: UInt,

    /// New-style CIGAR list (only used by splice-graph path, kept for parity).
    pub cigar: Vec<[u32; 2]>,

    pub var_ind: Vec<i32>,
    pub var_gen_coord: Vec<i32>,
    pub var_read_coord: Vec<i32>,
    pub var_allele: Vec<i8>,

    pub align_genes: std::collections::BTreeSet<u32>,
}

impl Default for Transcript {
    fn default() -> Self {
        Self::new()
    }
}

impl Transcript {
    /// C++ default ctor: zero everything.
    pub fn new() -> Self {
        Self {
            exons: vec![[0; EX_SIZE]; MAX_N_EXONS],
            shift_sj: vec![[0; 2]; MAX_N_EXONS],
            canon_sj: vec![0; MAX_N_EXONS],
            sj_annot: vec![0; MAX_N_EXONS],
            sj_str: vec![0; MAX_N_EXONS],
            intron_motifs: [0; 3],
            sj_motif_strand: 0,
            sj_yes: false,
            n_exons: 0,
            l_read: 0,
            read_length_pair_original: 0,
            i_read: 0,
            read_n_mates: 0,
            i_frag: 0,
            r_start: 0,
            ro_start: 0,
            r_length: 0,
            g_start: 0,
            g_length: 0,
            c_start: 0,
            chr: 0,
            str_: 0,
            ro_str: 0,
            haplo_type: 0,
            primary_flag: false,
            n_match: 0,
            n_mm: 0,
            mapped_length: 0,
            extend_l: 0,
            max_score: 0,
            n_gap: 0,
            l_gap: 0,
            n_del: 0,
            n_ins: 0,
            l_del: 0,
            l_ins: 0,
            n_unique: 0,
            n_anchor: 0,
            cigar: Vec::new(),
            var_ind: Vec::new(),
            var_gen_coord: Vec::new(),
            var_read_coord: Vec::new(),
            var_allele: Vec::new(),
            align_genes: std::collections::BTreeSet::new(),
        }
    }

    /// 1:1 port of `Transcript::reset()` (Transcript.cpp:8-26).
    ///
    /// Note: C++ only resets the scalar scoring/gap counters — the exon
    /// arrays (`exons`, `shiftSJ`, `canonSJ`, `sjAnnot`, `sjStr`) and
    /// `nExons` / `lRead` are preserved.
    pub fn reset(&mut self) {
        self.extend_l = 0;
        self.primary_flag = false;
        self.r_start = 0;
        self.ro_start = 0;
        self.r_length = 0;
        self.g_start = 0;
        self.g_length = 0;
        self.max_score = 0;
        self.n_match = 0;
        self.n_mm = 0;
        self.n_gap = 0;
        self.l_gap = 0;
        self.l_del = 0;
        self.l_ins = 0;
        self.n_del = 0;
        self.n_ins = 0;
        self.n_unique = 0;
        self.n_anchor = 0;
    }

    /// 1:1 port of `Transcript::add` (Transcript.cpp:28-36).
    pub fn add(&mut self, other: &Transcript) {
        self.max_score += other.max_score;
        self.n_match += other.n_match;
        self.n_mm += other.n_mm;
        self.n_gap += other.n_gap;
        self.l_gap += other.l_gap;
        self.l_del += other.l_del;
        self.n_del += other.n_del;
        self.l_ins += other.l_ins;
        self.n_ins += other.n_ins;
        self.n_unique += other.n_unique;
    }

    /// 1:1 port of `Transcript::extractSpliceJunctions` (Transcript.cpp:38-51).
    pub fn extract_splice_junctions(&self, sj_out: &mut Vec<[u64; 2]>) -> bool {
        let mut annot_yes = true;
        for iex in 0..self.n_exons.saturating_sub(1) as usize {
            if self.canon_sj[iex] >= 0 {
                let start = self.exons[iex][EX_G] + self.exons[iex][EX_L];
                let end = self.exons[iex + 1][EX_G] - start;
                sj_out.push([start, end]);
                if self.sj_annot[iex] == 0 {
                    annot_yes = false;
                }
            }
        }
        annot_yes
    }

    /// 1:1 port of `Transcript::chrStartLengthExtended` (Transcript.cpp:53-59).
    pub fn chr_start_length_extended(&self) -> u64 {
        let start1 = self.c_start - self.exons[0][EX_R];
        let last = self.n_exons as usize - 1;
        let length1 =
            self.exons[last][EX_G] + self.l_read - self.exons[last][EX_R] - self.exons[0][EX_G]
                + self.exons[0][EX_R];
        (start1 << 32) | length1
    }

    /// 1:1 port of `Transcript::alignScore` (Transcript_alignScore.cpp).
    ///
    /// Recomputes `max_score`, `n_mm`, `n_match` by re-scanning the exon
    /// sequence against `read1[0]` (forward) or `read1[2]` (RC) and
    /// re-adding junction / indel / gap scores. Matches the C++ loop
    /// exactly, including the `log2` genomic-length adjustment.
    pub fn align_score(
        &mut self,
        read1: &[Vec<u8>],
        g: &[u8],
        sjdb_score: i32,
        score_ins_base: i32,
        score_ins_open: i32,
        score_del_base: i32,
        score_del_open: i32,
        score_gap: i32,
        score_gap_gcag: i32,
        score_gap_atac: i32,
        score_gap_noncan: i32,
        score_genomic_length_log2_scale: f64,
    ) -> IntScore {
        self.max_score = 0;
        self.n_mm = 0;
        self.n_match = 0;
        if self.n_exons == 0 {
            return self.max_score;
        }
        let r: &[u8] = if self.ro_str == 0 {
            &read1[0]
        } else {
            &read1[2]
        };
        for iex in 0..self.n_exons as usize {
            let e = &self.exons[iex];
            for ii in 0..e[EX_L] {
                let r1 = r[(ii + e[EX_R]) as usize];
                let g1 = g[(ii + e[EX_G]) as usize];
                if r1 > 3 || g1 > 3 {
                    // skip
                } else if r1 == g1 {
                    self.max_score += 1;
                    self.n_match += 1;
                } else {
                    self.n_mm += 1;
                    self.max_score -= 1;
                }
            }
        }
        for iex in 0..self.n_exons.saturating_sub(1) as usize {
            if self.sj_annot[iex] == 1 {
                self.max_score += sjdb_score;
            } else {
                match self.canon_sj[iex] {
                    -3 => {}
                    -2 => {
                        let gap = self.exons[iex + 1][EX_R]
                            - self.exons[iex][EX_R]
                            - self.exons[iex][EX_L];
                        self.max_score += gap as i32 * score_ins_base + score_ins_open;
                    }
                    -1 => {
                        let gap = self.exons[iex + 1][EX_G]
                            - self.exons[iex][EX_G]
                            - self.exons[iex][EX_L];
                        self.max_score += gap as i32 * score_del_base + score_del_open;
                    }
                    0 => self.max_score += score_gap_noncan + score_gap,
                    1 | 2 => self.max_score += score_gap,
                    3 | 4 => self.max_score += score_gap_gcag + score_gap,
                    5 | 6 => self.max_score += score_gap_atac + score_gap,
                    _ => {}
                }
            }
        }
        if score_genomic_length_log2_scale != 0.0 {
            let last = self.n_exons as usize - 1;
            let span =
                (self.exons[last][EX_G] + self.exons[last][EX_L] - self.exons[0][EX_G]).max(1);
            let adj = ((span as f64).log2() * score_genomic_length_log2_scale - 0.5).ceil();
            self.max_score += adj as i32;
        }
        self.max_score
    }
}
