//! Ports `OutSJ.{h,cpp}` + `outputSJ.cpp` — per-chunk & global SJ collector
//! that writes `SJ.out.tab`.
//!
//! The C++ code uses a packed byte buffer (`Junction::dataSize` = 23). The
//! Rust port stores typed records and serializes them when needed; this keeps
//! the layout decoupled from memory alignment concerns.

use std::io::Write;

/// One splice junction recorded during alignment.
#[derive(Debug, Clone, Copy)]
pub struct Junction {
    pub start: u64,
    pub gap: u32,
    pub strand: u8,
    pub motif: u8,
    pub annot: u8,
    pub count_unique: u32,
    pub count_multiple: u32,
    pub overhang_left: u16,
    pub overhang_right: u16,
}

impl Junction {
    pub const DATA_SIZE: usize = 23;
}

/// 1:1 behavioral port of `class OutSJ` — a chunk-local collector for
/// junctions recorded during alignment.
#[derive(Debug, Default)]
pub struct OutSJ {
    pub junctions: Vec<Junction>,
}

impl OutSJ {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn n(&self) -> u64 {
        self.junctions.len() as u64
    }

    /// Port of `outputTranscriptSJ` (ReadAlign_outputTranscriptSJ.cpp:1-56).
    ///
    /// `tr_exons[iex] = (EX_G, EX_L)`; `tr_canon_sj`, `tr_sj_annot` are slices
    /// of the per-transcript arrays of length `n_exons`.
    pub fn record_transcript_sjs(
        &mut self,
        tr_exons: &[(u64, u64)],
        tr_canon_sj: &[i32],
        tr_sj_annot: &[u8],
        n_tr_out: u64,
        sj_read_start_n: u64,
    ) {
        if tr_exons.len() < 2 {
            return;
        }
        for iex in 0..tr_exons.len() - 1 {
            if tr_canon_sj[iex] < 0 {
                continue;
            }
            let start = tr_exons[iex].0 + tr_exons[iex].1;
            let gap = (tr_exons[iex + 1].0 - start) as u32;
            let overhang = std::cmp::min(tr_exons[iex].1, tr_exons[iex + 1].1) as u16;

            // De-duplicate per read (see C++ note about mate-overlap).
            let mut dup = false;
            for ii in sj_read_start_n as usize..self.junctions.len() {
                if self.junctions[ii].start == start && self.junctions[ii].gap == gap {
                    if self.junctions[ii].overhang_left < overhang {
                        self.junctions[ii].overhang_left = overhang;
                        self.junctions[ii].overhang_right = overhang;
                    }
                    dup = true;
                    break;
                }
            }
            if dup {
                continue;
            }

            let motif = tr_canon_sj[iex] as u8;
            let strand = if tr_canon_sj[iex] == 0 {
                0u8
            } else {
                ((tr_canon_sj[iex] + 1) % 2 + 1) as u8
            };
            let annot = tr_sj_annot[iex];
            let (count_unique, count_multiple) = if n_tr_out == 1 { (1, 0) } else { (0, 1) };

            self.junctions.push(Junction {
                start,
                gap,
                strand,
                motif,
                annot,
                count_unique,
                count_multiple,
                overhang_left: overhang,
                overhang_right: overhang,
            });
        }
    }

    /// Append `other` into `self`, in-order. Used after per-chunk chunks merge
    /// into a global `OutSJ` before the final `collapseSJ`.
    pub fn append(&mut self, other: &OutSJ) {
        self.junctions.extend_from_slice(&other.junctions);
    }

    /// Port of `Junction::outputStream` — writes one SJ.out.tab line
    /// (1-based, tab-separated). Caller must subtract chrStart from `start`
    /// and pick the correct chromosome name.
    pub fn write_junction_line(
        out: &mut impl Write,
        chr_name: &str,
        chr_start: u64,
        j: &Junction,
    ) -> std::io::Result<()> {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chr_name,
            j.start + 1 - chr_start,
            j.start + j.gap as u64 - chr_start,
            j.strand,
            j.motif,
            j.annot,
            j.count_unique,
            j.count_multiple,
            j.overhang_left,
        )
    }

    /// Port of `collapseSJ` — sort junctions and merge duplicates in place.
    ///
    /// C++ does this incrementally in a chunked buffer; we implement the same
    /// logic with a final sort-and-merge pass. Junctions sort by `(start,
    /// gap)` ascending. Merging sums `count_unique`/`count_multiple` and
    /// takes the max of overhangs.
    pub fn collapse_sj(&mut self) {
        self.junctions.sort_by(|a, b| {
            a.start
                .cmp(&b.start)
                .then_with(|| a.gap.cmp(&b.gap))
                .then_with(|| a.strand.cmp(&b.strand))
        });
        let mut collapsed: Vec<Junction> = Vec::with_capacity(self.junctions.len());
        for j in self.junctions.drain(..) {
            if let Some(last) = collapsed.last_mut() {
                if last.start == j.start && last.gap == j.gap && last.strand == j.strand {
                    last.count_unique += j.count_unique;
                    last.count_multiple += j.count_multiple;
                    last.overhang_left = last.overhang_left.max(j.overhang_left);
                    last.overhang_right = last.overhang_right.max(j.overhang_right);
                    if j.annot == 1 {
                        last.annot = 1;
                    }
                    continue;
                }
            }
            collapsed.push(j);
        }
        self.junctions = collapsed;
    }
}
