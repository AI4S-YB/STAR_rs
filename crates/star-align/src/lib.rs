//! star-align: read alignment core (M3/M4).
//!
//! Ports the following C++ sources (see comments per module):
//! - `Transcript.{h,cpp}`  -> [`transcript`]
//! - `ReadAlign.{h,cpp}`   -> [`read_align`]
//! - `ReadAlignChunk.{h,cpp}` + `ReadAlignChunk_processChunks.cpp` + `ReadAlignChunk_mapChunk.cpp`
//!    -> [`chunk`]
//! - `mapThreadsSpawn.cpp` -> [`threads`]
//! - `extendAlign.cpp`, `stitchAlignToTranscript.cpp`, `stitchWindowAligns.cpp`,
//!   `stitchGapIndel.cpp`, `binarySearch2.cpp`, `blocksOverlap.cpp`
//!   -> [`stitch`]
//! - `ClipMate_clip.cpp`, `ClipCR4.cpp` -> [`clip`]

pub mod calc_cigar;
pub mod chunk;
pub mod clip;
pub mod map_one_read;
pub mod mapped_filter;
pub mod mult_map_select;
pub mod one_read;
pub mod output_sam;
pub mod pe_overlap;
pub mod read_align;
pub mod rng;
pub mod seed;
pub mod stitch;
pub mod stitch_pieces;
pub mod stitch_window;
pub mod store_aligns;
pub mod threads;
pub mod transcript;
pub mod windows;

pub use read_align::ReadAlign;
pub use seed::{max_mappable_length_2strands, SeedHit};
pub use stitch_window::{stitch_window_aligns, StitchCtx};
pub use transcript::Transcript;
pub use windows::sj_align_split;
