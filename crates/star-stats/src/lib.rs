//! star-stats: ports `Stats.{h,cpp}` and the per-chunk/global stats totals
//! that back `Log.final.out`.

pub mod stats;

pub use stats::{Stats, TranscriptStatView, SJ_MOTIF_SIZE};
