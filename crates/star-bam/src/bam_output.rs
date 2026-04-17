//! Ports `BAMoutput.{h,cpp}` — per-thread BAM accumulator that bins records
//! by chromosome, flushes on threshold, and then concatenates per-bin files
//! into the final sorted BAM in M7.

// TODO(M7): pub struct BamOutput { ... }
