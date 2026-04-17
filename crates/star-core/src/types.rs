//! IncludeDefine.h: type aliases, exit codes, SAM / BAM / SJ / piece-coord constants.
//!
//! The C++ code defines:
//! ```cpp
//! #define uint  unsigned long long     // 64-bit
//! #define sint  signed long long
//! #define uint64 unsigned long long
//! #define uint32 unsigned int
//! #define uint16 unsigned short int
//! #define uchar  unsigned char
//! #define int64  long long
//! #define int32  int
//! #define uint128 __uint128_t
//! #define intScore int
//! ```
//! We use Rust primitive types but keep the same names via aliases for 1:1 mapping.

pub type UInt = u64;
pub type SInt = i64;
pub type UInt64 = u64;
pub type UInt32 = u32;
pub type UInt16 = u16;
pub type UChar = u8;
pub type Int64 = i64;
pub type Int32 = i32;
pub type UInt128 = u128;

pub type IntScore = i32;
pub type IntSwScore = i32;

pub const SCORE_MATCH: IntScore = 1;

// IncludeDefine.h
pub const GENOME_SPACING_CHAR: u8 = 5;

// Max 16-bit value used as a sentinel for window bin ("uintWinBinMax")
pub const UINT_WIN_BIN_MAX: u16 = u16::MAX;

// Compile-time feature selection: STARlong (COMPILE_FOR_LONG_READS) bumps these.
// Controlled via the `long-reads` Cargo feature; mirror the C++ `#define`.
#[cfg(not(feature = "long-reads"))]
pub const MAX_N_EXONS: usize = 20;
#[cfg(feature = "long-reads")]
pub const MAX_N_EXONS: usize = 1000;

#[cfg(not(feature = "long-reads"))]
pub const DEF_READ_SEQ_LENGTH_MAX: usize = 650;
#[cfg(feature = "long-reads")]
pub const DEF_READ_SEQ_LENGTH_MAX: usize = 500_000;

pub const MAX_N_MATES: usize = 3;
pub const DEF_READ_NAME_LENGTH_MAX: usize = 50_000;

pub const fn def_read_name_seq_length_max() -> usize {
    if DEF_READ_NAME_LENGTH_MAX > DEF_READ_SEQ_LENGTH_MAX {
        DEF_READ_NAME_LENGTH_MAX
    } else {
        DEF_READ_SEQ_LENGTH_MAX
    }
}

pub const BAM_OUTPUT_ONE_ALIGN_MAX_BYTES: usize = 100_000;

// SAM attribute ids (IncludeDefine.h)
pub const ATTR_NH: u16 = 1;
pub const ATTR_HI: u16 = 2;
pub const ATTR_AS: u16 = 3;
pub const ATTR_NM: u16 = 4;
pub const ATTR_MD: u16 = 5;
pub const ATTR_NM_LOWER: u16 = 6; // nM
pub const ATTR_JM: u16 = 7;
pub const ATTR_JI: u16 = 8;
pub const ATTR_XS: u16 = 9;
pub const ATTR_RG: u16 = 10;
pub const ATTR_VG: u16 = 11;
pub const ATTR_VA: u16 = 12;
pub const ATTR_VW: u16 = 13;
pub const ATTR_CH: u16 = 14;
pub const ATTR_MC: u16 = 15;
pub const ATTR_RB: u16 = 16;
pub const ATTR_CR: u16 = 17;
pub const ATTR_CY: u16 = 18;
pub const ATTR_UR: u16 = 19;
pub const ATTR_UY: u16 = 20;
pub const ATTR_CB: u16 = 21;
pub const ATTR_UB: u16 = 22;
pub const ATTR_GX: u16 = 23;
pub const ATTR_GN: u16 = 24;
pub const ATTR_SM: u16 = 25;
pub const ATTR_SS: u16 = 26;
pub const ATTR_SQ: u16 = 27;
pub const ATTR_HA: u16 = 28;
pub const ATTR_CN: u16 = 29;
pub const ATTR_GX_LOWER: u16 = 30;
pub const ATTR_GN_LOWER: u16 = 31;
pub const ATTR_SF: u16 = 32;

// BAM CIGAR op codes
pub const BAM_CIGAR_MAX_SIZE: usize = 10_000;
pub const BAM_CIGAR_OPERATION_SHIFT: u32 = 4;
pub const BAM_CIGAR_M: u8 = 0;
pub const BAM_CIGAR_I: u8 = 1;
pub const BAM_CIGAR_D: u8 = 2;
pub const BAM_CIGAR_N: u8 = 3;
pub const BAM_CIGAR_S: u8 = 4;
pub const BAM_CIGAR_H: u8 = 5;
pub const BAM_CIGAR_P: u8 = 6;
pub const BAM_CIGAR_EQ: u8 = 7;
pub const BAM_CIGAR_X: u8 = 8;

pub const BAM_ATTR_MAX_SIZE: usize = 10_000;

// Exit codes
pub const EXIT_CODE_BUG: i32 = 101;
pub const EXIT_CODE_PARAMETER: i32 = 102;
pub const EXIT_CODE_RUNTIME: i32 = 103;
pub const EXIT_CODE_INPUT_FILES: i32 = 104;
pub const EXIT_CODE_GENOME_FILES: i32 = 105;
pub const EXIT_CODE_SHM: i32 = 106;
pub const EXIT_CODE_GENOME_LOADING_WAITED_TOO_LONG: i32 = 107;
pub const EXIT_CODE_MEMORY_ALLOCATION: i32 = 108;
pub const EXIT_CODE_FILE_OPEN: i32 = 109;
pub const EXIT_CODE_FILE_WRITE: i32 = 110;
pub const EXIT_CODE_INCONSISTENT_DATA: i32 = 111;
pub const EXIT_CODE_FIFO: i32 = 112;

pub const EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS: i32 = 101;

pub const SJ_MOTIF_SIZE: usize = 7;
pub const SJ_SAM_ANNOTATED_MOTIF_SHIFT: u32 = 20;

pub const EXTEND_ORDER: u32 = 1;

pub const MAX_N_FRAG: usize = 2;
pub const MARK_FRAG_SPACER_BASE: u8 = 11;
pub const MAX_N_CHIMERAS: usize = 5;
pub const MAX_N_MULTMAP: usize = 100_000;
pub const MAX_SJ_REPEAT_SEARCH: usize = 255;
pub const MAX_QS_VALUE: u8 = 60;
pub const MAX_OUTPUT_FLAG: u32 = 10;

// Piece coord (PC_*) indices into a uiPC = [UInt; PC_SIZE]
pub const PC_R_START: usize = 0;
pub const PC_LENGTH: usize = 1;
pub const PC_STR: usize = 2;
pub const PC_DIR: usize = 3;
pub const PC_NREP: usize = 4;
pub const PC_SA_START: usize = 5;
pub const PC_SA_END: usize = 6;
pub const PC_I_FRAG: usize = 7;
pub const PC_SIZE: usize = 8;

pub const WC_STR: usize = 0;
pub const WC_CHR: usize = 1;
pub const WC_G_START: usize = 2;
pub const WC_G_END: usize = 3;
pub const WC_SIZE: usize = 4;

pub const WA_LENGTH: usize = 0;
pub const WA_R_START: usize = 1;
pub const WA_G_START: usize = 2;
pub const WA_NREP: usize = 3;
pub const WA_ANCHOR: usize = 4;
pub const WA_I_FRAG: usize = 5;
pub const WA_SJ_A: usize = 6;
pub const WA_SIZE: usize = 7;

pub const EX_R: usize = 0;
pub const EX_G: usize = 1;
pub const EX_L: usize = 2;
pub const EX_I_FRAG: usize = 3;
pub const EX_SJ_A: usize = 4;
pub const EX_SIZE: usize = 5;

// mapType
pub const MT_PE: usize = 0;
pub const MT_SIZE: usize = 5;

// Marker sentinel values in transcripts / alignments
pub const MARKER_ALL_PIECES_EXCEED_SEED_MULTIMAP_NMAX: u64 = 999_901;
pub const MARKER_NO_UNIQUE_PIECES: u64 = 999_902;
pub const MARKER_NO_GOOD_WINDOW: u64 = 999_903;
pub const MARKER_NO_GOOD_PIECES: u64 = 999_904;
pub const MARKER_TOO_MANY_ANCHORS_PER_WINDOW: u64 = 999_905;
pub const MARKER_MAX_N_MULT_EXCEEDED: u64 = 999_906;
pub const MARKER_FULL_LENGTH_MULTIMAPPER_EXCEEDED_ALIGN_WINDOWS_PER_READ_NMAX: u64 = 999_907;
pub const MARKER_ALL_PIECES_EXCEEDED_WIN_ANCHOR_MULTIMAP_NMAX: u64 = 999_908;
pub const MARKER_TOO_MANY_CHIMERAS: u64 = 999_909;
pub const MARKER_READ_TOO_SHORT: u64 = 999_910;

pub const PEMARKER_SINGLE_END: u32 = 0;
pub const PEMARKER_PAIR: u32 = 1;
pub const PEMARKER_ONE_END: u32 = 3;
pub const PEMARKER_TOO_MANY_PAIRS: u32 = 5;
pub const PEMARKER_CHIMERIC_PAIRS: u32 = 7;
pub const PEMARKER_CHIMERIC_SJ_READ1: u32 = 221;
pub const PEMARKER_CHIMERIC_SJ_READ2: u32 = 223;
pub const PEMARKER_CHIMERIC_SJ_READ1_AND_2: u32 = 225;
pub const PEMARKER_SINGLE_END_NOTMAPPED: u32 = 1001;

/// Pieces coords type (`uiPC`).
pub type UiPc = [UInt; PC_SIZE];
/// Window coords (`uiWC`).
pub type UiWc = [UInt; WC_SIZE];
/// Windowed alignment (`uiWA`).
pub type UiWa = [UInt; WA_SIZE];

/// STAR version literal; kept in one place.
///
/// Ports the `STAR_VERSION` macro from `source/VERSION`.
pub const STAR_VERSION: &str = "2.7.11b";
