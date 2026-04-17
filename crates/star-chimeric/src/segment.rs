//! 1:1 port of `ChimericSegment.{h,cpp}`.
//!
//! ChimericSegment is a thin value-type holding the read-coordinates
//! (roS / roE) and strand (`str_`) derived from a `Transcript` alignment
//! plus the per-mate read length (`readLength[0]`).
//!
//! In C++ `ChimericSegment` holds references to `Parameters` / `pCh` and
//! the `Transcript`. In Rust we decouple: `ChimericSegment::new` computes
//! the state, and [`segment_check`] evaluates the quality gate against the
//! `ParametersChimeric` threshold.

use star_align::transcript::Transcript;
use star_core::{EX_L, EX_R, UInt};
use star_params::parameters::ParametersChimeric;

/// Port of `class ChimericSegment`.
///
/// - `str_`: 0 = undefined; 1 = same as RNA; 2 = opposite to RNA.
/// - `ro_s` / `ro_e`: start/end in the read-original coordinate system
///   (covers the read sequence as read from the FASTQ, regardless of the
///   alignment strand).
#[derive(Clone, Debug)]
pub struct ChimericSegment {
    pub str_: u64,
    pub ro_s: UInt,
    pub ro_e: UInt,
}

impl ChimericSegment {
    /// Port of the `ChimericSegment::ChimericSegment(...)` ctor.
    ///
    /// `read_length_0` is the per-mate read length (`readLength[0]`). In SE
    /// and in PE before mate-split it corresponds to the length of mate 1.
    pub fn new(align: &Transcript, read_length_0: UInt) -> Self {
        let im = align.intron_motifs;
        let str_ = if (im[1] == 0 && im[2] == 0) || (im[1] > 0 && im[2] > 0) {
            0
        } else if (align.str_ == 0) == (im[1] > 0) {
            1
        } else {
            2
        };

        let n_exons = align.n_exons as usize;
        let last = n_exons.saturating_sub(1);
        let l_read = align.l_read;

        let mut ro_s = if align.str_ == 0 {
            align.exons[0][EX_R]
        } else {
            l_read - align.exons[last][EX_R] - align.exons[last][EX_L]
        };
        let mut ro_e = if align.str_ == 0 {
            align.exons[last][EX_R] + align.exons[last][EX_L] - 1
        } else {
            l_read - align.exons[0][EX_R] - 1
        };
        if ro_s > read_length_0 {
            ro_s -= 1;
        }
        if ro_e > read_length_0 {
            ro_e -= 1;
        }

        Self { str_, ro_s, ro_e }
    }
}

/// Port of `ChimericSegment::segmentCheck`.
pub fn segment_check(align: &Transcript, p_ch: &ParametersChimeric) -> bool {
    align.r_length >= p_ch.segment_min && align.intron_motifs[0] == 0
}
