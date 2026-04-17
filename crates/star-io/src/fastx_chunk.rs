//! Port of `ReadAlignChunk_processChunks_readChunkFastx.cpp`.
//!
//! Splits streamed FASTX input into per-thread byte chunks bounded by
//! `--readMapNumber`/`--runThreadN`. For M3 single-threaded we read
//! whole-file into a `Vec<u8>` and offer a simple iterator; M4 wires this
//! into per-thread buffers with size-matched chunk grain.

// TODO(M3): pub fn read_chunk_fastx(...) -> usize (bytes copied).
