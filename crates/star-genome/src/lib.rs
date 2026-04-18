//! star-genome: ports `Genome.{h,cpp}`, `Genome_genomeLoad.cpp`,
//! `Genome_genomeGenerate.cpp`, `genomeSAindex.cpp`, `genomeScanFastaFiles.cpp`,
//! `SuffixArrayFuns.cpp`, `genomeParametersWrite.cpp`, and sjdb helpers.
//!
//! Fully ported in M2 (real code):
//! - [`genome::Genome`] struct wires the core state (G buffer, chrStart/Length,
//!   chrName, nChrReal, PackedArray handles, sjdb fields, ...).
//! - [`fasta::genome_scan_fasta_files`] reproduces `genomeScanFastaFiles`.
//! - [`io::write_chr_info`] reproduces `Genome::writeChrInfo`.
//! - [`io::write_genome_sequence`] reproduces `Genome::writeGenomeSequence`.
//! - [`io::genome_parameters_write`] reproduces `genomeParametersWrite.cpp`.
//! - [`sa::compare_suffixes`] reproduces the `funCompareSuffixes` qsort cmp.
//!
//! Still stubbed (clearly marked): SA sort driver, SAindex, sjdb.

pub mod fasta;
pub mod generate;
pub mod genome;
pub mod io;
pub mod load;
pub mod sa;
pub mod sa_funs;
pub mod sa_index;
pub mod sa_sort;

pub use generate::genome_generate;
pub use genome::Genome;
pub use load::genome_load;
