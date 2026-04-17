//! star-quant: transcriptome BAM + gene counts.
//!
//! Ports (module → C++ source):
//! - [`transcriptome`] -> `Transcriptome.{h,cpp}`
//!                        + `Transcriptome_geneFullAlignOverlap.cpp`
//!                        + `Transcriptome_quantAlign.cpp`
//! - [`quant`]         -> `Quantifications.{h,cpp}`
//!                        + `ReadAlign_quantTranscriptome.cpp`
//!                        + `ReadAlign_outputTranscriptCIGAR.cpp` (transcriptome)
//!
//! M7 populates these modules.

pub mod quant;
pub mod transcriptome;
