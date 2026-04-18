//! star-bam: SAM/BAM output for STAR.
//!
//! Ports (module → C++ source):
//! - [`sam_headers`]       -> `samHeaders.cpp`
//! - [`sam_output`]        -> `outputSJ.cpp` + `ReadAlign_outputTranscriptSAM.cpp` helpers
//! - [`bam_output`]        -> `BAMoutput.{h,cpp}`
//! - [`bam_functions`]     -> `BAMfunctions.{h,cpp}`
//! - [`bam_sort_coord`]    -> `bamSortByCoordinate.cpp` + `BAMbinSortByCoordinate.cpp`
//!   + `BAMbinSortUnmapped.cpp`
//! - [`bam_cat`]           -> `bam_cat.c`
//! - [`bam_remove_dups`]   -> `bamRemoveDuplicates.cpp`
//! - [`signal_from_bam`]   -> `signalFromBAM.cpp`
//!
//! M3 fills in `sam_output` + `sam_headers` only. M7 completes BAM paths.

pub mod bam_cat;
pub mod bam_functions;
pub mod bam_output;
pub mod bam_remove_dups;
pub mod bam_sort_bin;
pub mod bam_sort_coord;
pub mod bam_sort_feed;
pub mod bam_sort_record;
pub mod bam_sort_unmapped;
pub mod sam_headers;
pub mod sam_output;
pub mod sam_to_bam;
pub mod signal_from_bam;
