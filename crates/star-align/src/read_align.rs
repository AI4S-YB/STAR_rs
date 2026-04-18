//! Ports `class ReadAlign` (ReadAlign.h) state and methods.
//!
//! The struct mirrors the C++ layout with types translated to idiomatic Rust
//! (`Vec` / `Option`). Public fields intentionally match C++ names (snake
//! case) so the ported method bodies can be cross-checked line-by-line.

use crate::rng::Mt19937;
use crate::transcript::Transcript;
use star_core::types::{IntScore, MAX_N_MATES, PC_SIZE, UINT_WIN_BIN_MAX, WA_SIZE, WC_SIZE};
use star_params::parameters::Parameters;

/// Piece-coordinate entry (`uiPC` in C++).
pub type PieceCoord = [u64; PC_SIZE];
/// Window-coordinate entry (`uiWC`).
pub type WinCoord = [u64; WC_SIZE];
/// Window-alignment entry (`uiWA`).
pub type WinAlign = [u64; WA_SIZE];

/// Port of `class ReadAlign` (ReadAlign.h:22-235).
///
/// References to `Parameters`, `Genome`, `Transcriptome`, `BAMoutput`,
/// `OutSJ`, `SoloRead`, `SpliceGraph` etc. are held externally by the call
/// sites; we pass them as borrowed arguments to individual methods. This
/// keeps `ReadAlign` a pure state struct, avoiding self-referential borrow
/// problems while matching the C++ control flow.
pub struct ReadAlign {
    // Book-keeping.
    pub i_read: u64,
    pub i_read_all: u64,
    pub read_files_index: u32,
    pub read_filter: u8,
    pub read_file_type: i32,
    pub revert_strand: bool,

    // PRNGs.
    pub rng: Mt19937,
    pub rng_mult_order: Mt19937,

    // Read data — `Read1[iMate]` is ASCII text; internal numeric form is
    // stored in `read0` and reverse-complement in `read2` (collapsed to
    // `read_buf[0..3]` for forward / padded-forward / reverse-complement).
    pub read_name: String,
    pub read_name_mates: Vec<String>,
    pub read_name_extra: Vec<String>,
    pub read1: Vec<Vec<u8>>, // 0/1/2: FWD ASCII / PAD ASCII / RC ASCII
    pub read0: Vec<Vec<u8>>, // numeric per mate
    pub qual0: Vec<Vec<u8>>, // quality per mate
    pub read_length: [u64; MAX_N_MATES],
    pub read_length_original: [u64; MAX_N_MATES],
    pub read_length_pair: u64,
    pub read_length_pair_original: u64,
    pub l_read: u64,

    // Per-mate max score snapshot.
    pub max_score_mate: [IntScore; MAX_N_MATES],
    pub max_score: IntScore,

    // Seeds / pieces.
    pub pc: Vec<PieceCoord>,
    pub n_p: u64,
    pub n_a: u64,
    pub n_um: [u64; 2],

    /// `splitR[3][maxNsplit]` in C++; stored column-major as `[start,length,frag]`.
    pub split_r: Vec<[u64; 3]>,
    pub n_split: u64,
    pub read_n_mates: u32,

    // Sort buffers.
    pub stored_lmin: u64,
    pub uniq_lmax: u64,
    pub uniq_lmax_ind: u64,
    pub mult_lmax: u64,
    pub mult_lmax_n: u64,
    pub mult_nmin: u64,
    pub mult_nmin_l: u64,
    pub mult_nmax: u64,
    pub mult_nmax_l: u64,

    // Windows.
    pub wc: Vec<WinCoord>,
    pub win_bin: [Vec<u16>; 2],
    pub n_w: u64,
    pub n_wall: u64,

    // Aligns per window.
    pub wa: Vec<Vec<WinAlign>>,
    pub n_wa: Vec<u64>,
    pub n_wap: Vec<u64>,
    pub wa_lrec: Vec<u64>,
    pub wlast_anchor: Vec<u64>,
    pub wa_incl: Vec<bool>,

    // Long-reads SW window coverage buffers.
    pub sw_win_cov: Vec<u64>,

    // Transcripts.
    pub tr_init: Transcript,
    pub tr_a: Transcript,
    pub tr_a1: Transcript,
    pub tr_array: Vec<Transcript>,
    /// `trAll[iW][iTr]`; Rust port stores transcripts inline rather than
    /// holding pointers into `tr_array`.
    pub tr_all: Vec<Vec<Transcript>>,
    pub n_win_tr: Vec<u64>,
    pub tr_best: u64, // window index; UINT64_MAX if unset
    pub tr_mult: Vec<u64>,
    pub n_tr: u64,

    // Misc.
    pub map_marker: u64,
    pub unmap_type: i32,
    pub out_filter_mismatch_nmax_total: u64,
    pub mate_mapped: [bool; 2],

    /// Indices into `tr_array` (when non-empty) OR inline storage — we keep
    /// it inline to avoid cross-reference lifetime issues.
    pub tr_mult_array: Vec<Transcript>,
    /// Per-mate clipping left/right (stubbed; lengths only).
    pub clip_mates: [[ClippedMate; 2]; 2],
    /// Buffered CIGAR strings per mate (produced by `calc_cigar`).
    pub mates_cigar: Vec<String>,

    /// Port of `ReadAlign::peOv` — paired-end mate-overlap state.
    pub pe_ov: PeOverlap,

    // --- M6 chimeric state (ReadAlign.h:153-168 fragment) ---
    /// `Transcript trChim[MAX_CHIMNUM]`; we only use index 0/1.
    pub tr_chim: [Transcript; 2],
    /// Set to true when a chimeric alignment is recorded.
    pub chim_record: bool,
    /// Number of chimeric segments actually recorded (0 or 2).
    pub chim_n: u64,
    /// Chimeric junction coordinates (genome, 0-based).
    pub chim_j0: u64,
    pub chim_j1: u64,
    /// Chimeric junction motif: 0 non-canonical, 1 GT-AG, 2 CT-AC, -1 mate-bracketed.
    pub chim_motif: i32,
    /// Chimeric strand: 0 undefined, 1 same-as-RNA, 2 opposite.
    pub chim_str: u64,
    /// Length of repeat flanking the chimeric junction.
    pub chim_repeat0: u64,
    pub chim_repeat1: u64,
}

/// Port of the anonymous struct inside `ReadAlign::peOv`.
#[derive(Default, Clone, Copy, Debug)]
pub struct PeOverlap {
    pub yes: bool,
    pub n_ov: u64,
    pub ov_s: u64,
    pub mate_start: [u64; 2],
}

/// Minimal port of `ClipMate` — only `clippedN` is tracked in M3.
#[derive(Default, Clone, Copy, Debug)]
pub struct ClippedMate {
    pub clipped_n: u64,
}

const UNSET_TR: u64 = u64::MAX;

impl ReadAlign {
    pub fn new(run_rng_seed: u32) -> Self {
        Self {
            i_read: 0,
            i_read_all: 0,
            read_files_index: 0,
            read_filter: b'N',
            read_file_type: 0,
            revert_strand: false,
            rng: Mt19937::new(run_rng_seed),
            rng_mult_order: Mt19937::new(run_rng_seed),
            read_name: String::new(),
            read_name_mates: Vec::new(),
            read_name_extra: Vec::new(),
            read1: Vec::new(),
            read0: Vec::new(),
            qual0: Vec::new(),
            read_length: [0; MAX_N_MATES],
            read_length_original: [0; MAX_N_MATES],
            read_length_pair: 0,
            read_length_pair_original: 0,
            l_read: 0,
            max_score_mate: [0; MAX_N_MATES],
            max_score: 0,
            pc: Vec::new(),
            n_p: 0,
            n_a: 0,
            n_um: [0, 0],
            split_r: Vec::new(),
            n_split: 0,
            read_n_mates: 1,
            stored_lmin: 0,
            uniq_lmax: 0,
            uniq_lmax_ind: 0,
            mult_lmax: 0,
            mult_lmax_n: 0,
            mult_nmin: 0,
            mult_nmin_l: 0,
            mult_nmax: 0,
            mult_nmax_l: 0,
            wc: Vec::new(),
            win_bin: [Vec::new(), Vec::new()],
            n_w: 0,
            n_wall: 0,
            wa: Vec::new(),
            n_wa: Vec::new(),
            n_wap: Vec::new(),
            wa_lrec: Vec::new(),
            wlast_anchor: Vec::new(),
            wa_incl: Vec::new(),
            sw_win_cov: Vec::new(),
            tr_init: Transcript::new(),
            tr_a: Transcript::new(),
            tr_a1: Transcript::new(),
            tr_array: Vec::new(),
            tr_all: Vec::new(),
            n_win_tr: Vec::new(),
            tr_best: UNSET_TR,
            tr_mult: Vec::new(),
            n_tr: 0,
            map_marker: 0,
            unmap_type: 0,
            out_filter_mismatch_nmax_total: 0,
            mate_mapped: [false; 2],
            tr_mult_array: Vec::new(),
            clip_mates: [[ClippedMate::default(); 2]; 2],
            mates_cigar: Vec::new(),
            pe_ov: PeOverlap::default(),
            tr_chim: [Transcript::new(), Transcript::new()],
            chim_record: false,
            chim_n: 0,
            chim_j0: 0,
            chim_j1: 0,
            chim_motif: 0,
            chim_str: 0,
            chim_repeat0: 0,
            chim_repeat1: 0,
        }
    }

    /// Allocate buffers whose size depends on `Parameters`. Mirrors the
    /// `new uintWinBin[P.winBinN]` / seed-buffer allocations in the C++
    /// `ReadAlign::ReadAlign` constructor.
    pub fn init_from_params(&mut self, p: &Parameters) {
        let wbn = p.win_bin_n as usize;
        self.win_bin[0] = vec![UINT_WIN_BIN_MAX; wbn];
        self.win_bin[1] = vec![UINT_WIN_BIN_MAX; wbn];

        let wpr = p.align_window_per_read_nmax as usize;
        let tpw = p.align_transcripts_per_window_nmax as usize;
        self.wc = vec![[0u64; WC_SIZE]; wpr];
        self.wa = vec![vec![[0u64; WA_SIZE]; p.seed_per_window_nmax as usize]; wpr];
        self.n_wa = vec![0u64; wpr];
        self.n_wap = vec![0u64; wpr];
        self.wa_lrec = vec![0u64; wpr];
        self.wlast_anchor = vec![0u64; wpr];
        self.wa_incl = vec![false; p.seed_per_window_nmax as usize];
        self.n_win_tr = vec![0u64; wpr];
        self.tr_all = vec![Vec::with_capacity(tpw); wpr];

        let pc_cap = p.seed_per_read_nmax as usize;
        self.pc = vec![[0u64; PC_SIZE]; pc_cap];
        self.split_r = vec![[0u64; 3]; p.max_n_split as usize];

        self.read_n_mates = p.read_nmates;
    }

    /// Port of `ReadAlign::resetN` (ReadAlign.cpp). Zero per-read counters
    /// keeping buffers allocated.
    pub fn reset_n(&mut self) {
        self.n_p = 0;
        self.n_a = 0;
        self.n_um = [0, 0];
        self.n_w = 0;
        self.n_wall = 0;
        self.n_tr = 0;
        self.map_marker = 0;
        self.unmap_type = -1;
        self.max_score = 0;
        self.max_score_mate = [0; MAX_N_MATES];
        self.stored_lmin = 0;
        self.uniq_lmax = 0;
        self.uniq_lmax_ind = 0;
        self.mult_lmax = 0;
        self.mult_lmax_n = 0;
        self.mult_nmin = 0;
        self.mult_nmin_l = 0;
        self.mult_nmax = 0;
        self.mult_nmax_l = 0;
        self.tr_best = UNSET_TR;
        self.chim_record = false;
        self.chim_n = 0;
        self.chim_j0 = 0;
        self.chim_j1 = 0;
        self.chim_motif = 0;
        self.chim_str = 0;
        self.chim_repeat0 = 0;
        self.chim_repeat1 = 0;
    }
}
