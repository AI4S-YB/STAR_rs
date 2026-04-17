//! Ports `BAMfunctions.{h,cpp}` — BAM primitive helpers (pack CIGAR, encode
//! sequence, compute bin, etc.). Backed by `noodles-bam` primitives where
//! possible, with thin wrappers matching the C++ signatures for 1:1 porting.

// TODO(M7): pub fn bam_reg2bin(...) / bam_pack_cigar(...) ...
