//! 1:1 port of `class Genome` state (see `Genome.h`).

use std::collections::HashMap;

use star_core::packed::PackedArray;
use star_core::types::GENOME_SPACING_CHAR;

/// Byte offset of `G` inside the owning `g1` allocation.
///
/// Matches C++ `L = 200` in `Genome_genomeLoad.cpp`: STAR reserves 200 bytes
/// of pre-padding and 200 bytes of post-padding so that seed-extension
/// readers (up to a mate length, typically ≤ 300 bytes) can safely cross
/// either end of the genome sequence without overrunning. The matching
/// pre/post-pad size is also used during `genomeGenerate` via
/// [`Genome::genome_sequence_allocate`] so the SA/SAi built against an
/// in-memory genome can be loaded back at the identical offset.
pub const G_OFFSET: usize = 200;

/// Mirrors `Genome.h` field layout (omitting Variation / SuperTranscriptome
/// for milestones ≤ M7).
#[derive(Clone, Default)]
pub struct Genome {
    /// Full genome sequence (numeric 0..=4) of length `n_genome`. In the C++
    /// code this is a `char*`; we own a `Vec<u8>` and expose slices.
    pub g: Vec<u8>,
    /// `char *G1`: allocation backing; when `--genomeLoad` uses shared memory,
    /// `G1 = shmStart` and `G = G1 + shmHeader`. For M2 we only need the
    /// in-process backing, so `g1_alloc` is the allocation length.
    pub n_g1_alloc: u64,
    pub n_genome: u64,

    pub sa: PackedArray,
    pub sa_insert: PackedArray,
    pub sa_pass1: PackedArray,
    pub sa_pass2: PackedArray,
    pub sai: PackedArray,

    pub n_genome_insert: u64,
    pub n_genome_pass1: u64,
    pub n_genome_pass2: u64,
    pub n_sa_insert: u64,
    pub n_sa_pass1: u64,
    pub n_sa_pass2: u64,

    // Chromosomes
    pub chr_start: Vec<u64>,
    pub chr_length: Vec<u64>,
    pub chr_length_all: Vec<u64>,
    pub chr_name: Vec<String>,
    pub chr_name_all: Vec<String>,
    pub chr_name_index: HashMap<String, u64>,

    pub genome_chr_bin_nbases: u64,
    pub chr_bin_n: u64,
    pub chr_bin: Vec<u64>,

    pub genome_sa_index_start: Vec<u64>,

    pub n_sa: u64,
    pub n_sa_byte: u64,
    pub n_chr_real: u64,

    pub n_sai: u64,
    pub g_strand_bit: u8,
    pub sai_mark_n_bit: u8,
    pub sai_mark_absent_bit: u8,
    pub g_strand_mask: u64,
    pub sai_mark_absent_mask: u64,
    pub sai_mark_absent_mask_c: u64,
    pub sai_mark_n_mask: u64,
    pub sai_mark_n_mask_c: u64,

    // SJ database state
    pub sjdb_overhang: u64,
    pub sjdb_length: u64,
    pub sj_chr_start: u64,
    pub sjdb_n: u64,
    pub sj_g_start: u64,
    pub sjdb_start: Vec<u64>,
    pub sjdb_end: Vec<u64>,
    pub sj_d_start: Vec<u64>,
    pub sj_a_start: Vec<u64>,
    pub sj_str: Vec<u64>,
    pub sjdb_motif: Vec<u8>,
    pub sjdb_shift_left: Vec<u8>,
    pub sjdb_shift_right: Vec<u8>,
    pub sjdb_strand: Vec<u8>,

    pub genome_insert_l: u64,
    pub genome_insert_chr_ind_first: u64,
}

impl Genome {
    pub fn new() -> Self {
        Self::default()
    }

    /// Port of `Genome::genomeSequenceAllocate`:
    ///
    /// ```cpp
    /// nG1allocOut = (nGenomeIn + 100) * 2;
    /// G1out = new char[nG1allocOut];
    /// Gout  = G1out + 100;
    /// memset(G1out, GENOME_spacingChar, nG1allocOut);
    /// ```
    ///
    /// We allocate the full `g1` buffer, pre-fill with `GENOME_spacingChar`
    /// (== 5), and set `n_genome`. Callers index through [`Self::g_slice`]
    /// or the `g` Vec directly using `G_OFFSET`.
    ///
    /// We allocate a bit more than C++ to keep the `G_OFFSET = 200` pre-pad
    /// consistent with the loader (see [`G_OFFSET`] docs).
    pub fn genome_sequence_allocate(&mut self, n_genome: u64) {
        self.n_genome = n_genome;
        self.n_g1_alloc = (n_genome + G_OFFSET as u64) * 2;
        self.g = vec![GENOME_SPACING_CHAR; self.n_g1_alloc as usize];
    }

    /// Immutable pointer to `G[0]` (i.e. past the 100-byte pre-padding).
    ///
    /// # Safety
    /// Returned pointer is valid for reads of the full G view
    /// `G[-100 .. n_g1_alloc - 100]` (so `G[-3]` in SA loops is fine).
    pub fn g_ptr(&self) -> *const u8 {
        unsafe { self.g.as_ptr().add(G_OFFSET) }
    }

    pub fn g_ptr_mut(&mut self) -> *mut u8 {
        unsafe { self.g.as_mut_ptr().add(G_OFFSET) }
    }

    /// `chrBinFill` port (see `Genome.cpp`): builds `chr_bin[0..chr_bin_n]`
    /// that maps each `genomeChrBinNbases`-sized bin to its chr index.
    pub fn chr_bin_fill(&mut self) {
        self.chr_bin_n = self.n_genome / self.genome_chr_bin_nbases + 1;
        self.chr_bin = vec![0u64; self.chr_bin_n as usize];
        for ii in 0..self.n_chr_real as usize {
            let b1 = (self.chr_start[ii] / self.genome_chr_bin_nbases) as usize;
            let b2 = (self.chr_start[ii + 1] / self.genome_chr_bin_nbases) as usize;
            for jj in b1..b2 {
                self.chr_bin[jj] = ii as u64;
            }
        }
    }
}
