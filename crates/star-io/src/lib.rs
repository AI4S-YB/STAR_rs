//! star-io: input/output streams & read loading.
//!
//! Ports:
//! - `InOutStreams.{h,cpp}` -> [`streams`]
//! - `readLoad.{h,cpp}`     -> [`read_load`]
//! - `Parameters_openReadsFiles.cpp` + `..._closeReadsFiles.cpp`
//!   + `..._readFilesInit.cpp` + `..._readSAMheader.cpp` -> [`reads_files`]
//! - `ReadAlignChunk_processChunks_readChunkFastx.cpp` -> [`fastx_chunk`]

pub mod fastx_chunk;
pub mod read_load;
pub mod reads_files;
pub mod sam_headers;
pub mod streams;
