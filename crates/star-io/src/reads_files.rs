//! Ports:
//! - `Parameters_openReadsFiles.cpp`
//! - `Parameters_closeReadsFiles.cpp`
//! - `Parameters_readFilesInit.cpp`
//! - `Parameters_readSAMheader.cpp`
//!
//! Opens `--readFilesIn` inputs (optionally through a `--readFilesCommand`
//! pipe) and wires them into `InOutStreams`. SAM inputs require extracting
//! the header first to write @SQ lines to the output SAM.

// TODO(M3): open_reads_files(params, io) -> anyhow::Result<()>
// TODO(M3): close_reads_files(io) -> ()
// TODO(M3): read_files_init(params) -> ()
// TODO(M3): read_sam_header(params, io) -> ()
