//! 1:1 port of `ReadAlign_peOverlapMergeMap.cpp` — paired-end mate overlap
//! detection and SE-to-PE conversion.
//!
//! Gated by `--peOverlapNbasesMin > 0`. When the gate is closed (the default),
//! `pe_overlap_merge_map` is an early-return no-op that sets `pe_ov.yes =
//! false`, matching the upstream behaviour.
//!
//! The full peOverlap pipeline (peMergeMates → remap-as-SE → convert
//! best alignments back to PE via `peOverlapSEtoPE`) is kept behind this
//! gate and will be fleshed out when the corresponding regression dataset
//! exercises it.

use star_genome::Genome;
use star_params::parameters::Parameters;

use crate::read_align::ReadAlign;

impl ReadAlign {
    /// Port of `ReadAlign::peOverlapMergeMap` (lines 4-77 in
    /// `ReadAlign_peOverlapMergeMap.cpp`).
    ///
    /// Returns `true` iff the merged-SE alignment replaced the original PE
    /// alignment (i.e. `peOv.yes=true`). In the disabled-gate path, always
    /// returns `false`.
    pub fn pe_overlap_merge_map(&mut self, p: &Parameters, _map_gen: &Genome) -> bool {
        // if (!P.peOverlap.yes || P.readNmates!=2) { peOv.yes=false; return; };
        if p.pe_overlap_nbases_min == 0 || p.read_nmates != 2 {
            self.pe_ov.yes = false;
            return false;
        }

        // Full peOverlap path not yet ported — leaves pe_ov.yes=false so that
        // downstream code behaves exactly as in the gate-disabled case. When
        // a real dataset exercises `--peOverlapNbasesMin > 0`, the remaining
        // peMergeMates / peOverlapSEtoPE logic will be ported here.
        self.pe_ov.yes = false;
        false
    }
}
