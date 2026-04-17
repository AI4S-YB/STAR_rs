//! Adapter clipping: ports `ClipMate_clip.cpp` and `ClipCR4.cpp`.
//!
//! `ClipCR4` delegates to opal via `opal-sys` FFI for AVX2 Smith-Waterman.
//! M3 provides entry points; real logic to be filled during alignment port.

// TODO(M3): pub struct ClipMate; impl ClipMate { pub fn clip(&mut self, ...) }
// TODO(M3): pub fn clip_cr4_search(...) -> calls opal_sys::opalSearchDatabase.
