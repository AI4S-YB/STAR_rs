//! star-sjdb: splice-junction database + two-pass support.
//!
//! Ports (module → C++ source):
//! - [`prepare`]     -> `sjdbPrepare.cpp` + `sjdbLoadFromFiles.cpp`
//!                      + `sjdbLoadFromStream.cpp` + `sjdbReadBED.cpp`
//!                      + `Genome_insertSequences.cpp`
//!                      + `sjdbInsertJunctions.cpp` + `sjdbBuildIndex.cpp`
//! - [`out_sj`]      -> `OutSJ.{h,cpp}` + `outputSJ.cpp`
//! - [`twopass`]     -> `twoPassRunPass1.cpp`
//! - [`loader`]      -> `sjAlignSplit.cpp`
//!
//! M5 populates these modules; M2 only emits a minimal `sjdb/` shell.

pub mod binary_search2;
pub mod build_index;
pub mod insert_junctions;
pub mod loader;
pub mod out_sj;
pub mod output_sj;
pub mod prepare;
pub mod sjdb_class;
pub mod twopass;

pub use sjdb_class::SjdbLoci;
