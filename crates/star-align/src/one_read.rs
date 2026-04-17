//! 1:1 port of `ReadAlign::oneRead` (ReadAlign_oneRead.cpp).
//!
//! Handles: load from buffered stream, combine mates (PE), build
//! complement/reverse buffers, stats, then call `map_one_read` and downstream
//! output.
//!
//! Output routing (BAM / chimeric / WASP) is delegated to callbacks the caller
//! provides, so this module stays independent of star-io / star-chimeric.

use star_core::seq::complement_seq_numbers;
use star_core::types::{DEF_READ_SEQ_LENGTH_MAX, MARK_FRAG_SPACER_BASE};
use star_genome::Genome;
use star_io::read_load::{read_load, ClipMate, LoadedRead};
use star_params::parameters::Parameters;
use star_stats::stats::Stats;

use crate::read_align::ReadAlign;

/// Paths `oneRead()` can take for output / chimeric detection / wasp; stubbed
/// for M3 basic SE path.
#[derive(Debug, Default)]
pub struct OneReadOutputs {
    pub mapped: bool,
    pub is_chimeric: bool,
}

/// Port of `ReadAlign::oneRead` return-code convention (ReadAlign.h:129-136).
///
/// - `ReadDone` (0): read was loaded and processed.
/// - `NoMoreReads` (-1): EOF — no read consumed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OneReadStatus {
    ReadDone,
    NoMoreReads,
}

impl ReadAlign {
    /// Drive one read through the pipeline. Returns `Ok(None)` on EOF.
    pub fn one_read<R: std::io::BufRead>(
        &mut self,
        p: &Parameters,
        map_gen: &Genome,
        readers: &mut [R],
        stats: &mut Stats,
    ) -> anyhow::Result<Option<OneReadOutputs>> {
        // Load each mate.
        let mut loaded: Vec<LoadedRead> = Vec::with_capacity(readers.len());
        for reader in readers.iter_mut() {
            let mut clip = [ClipMate::default(), ClipMate::default()];
            match read_load(reader, p, &mut clip)? {
                Some(rec) => loaded.push(rec),
                None => return Ok(None),
            }
        }
        self.one_read_loaded(p, map_gen, loaded, stats)
    }

    /// Variant of `one_read` that accepts already-loaded mate records. Used
    /// by multi-threaded drivers that pre-fetch records on a serial reader
    /// thread and dispatch them to workers.
    pub fn one_read_loaded(
        &mut self,
        p: &Parameters,
        map_gen: &Genome,
        loaded: Vec<LoadedRead>,
        stats: &mut Stats,
    ) -> anyhow::Result<Option<OneReadOutputs>> {
        if loaded.is_empty() {
            return Ok(None);
        }
        // Drop read name / quality into ReadAlign state.
        self.read_name = loaded[0].read_name.clone();
        self.read_name_mates = loaded.iter().map(|r| r.read_name.clone()).collect();
        self.read_name_extra = loaded.iter().map(|r| r.read_name_extra.clone()).collect();
        self.read_files_index = loaded[0].read_files_index;
        self.i_read_all = loaded[0].i_read_all;
        self.read_filter = loaded[0].read_filter as u8;
        self.read_file_type = loaded[0].file_type;
        for (i, rec) in loaded.iter().enumerate() {
            if i < self.read_length.len() {
                self.read_length[i] = rec.lread;
                self.read_length_original[i] = rec.lread_original;
            }
        }

        // Ensure Read0 / Read1 / Qual0 are allocated with capacity 3 entries.
        while self.read0.len() < 3 {
            self.read0.push(Vec::new());
        }
        while self.read1.len() < 3 {
            self.read1.push(Vec::new());
        }
        while self.qual0.len() < 3 {
            self.qual0.push(Vec::new());
        }

        // Populate Read0 / Qual0 with the ASCII sequence / quality per mate
        // (matches C++ `Read0[iMate]` and `Qual0[iMate]`).
        for (i, rec) in loaded.iter().enumerate() {
            if i < self.read0.len() {
                self.read0[i].clear();
                self.read0[i].extend_from_slice(&rec.seq);
                self.qual0[i].clear();
                self.qual0[i].extend_from_slice(&rec.qual);
            }
        }

        // Combine mates into a single Lread buffer in `read1[0]` + `read1[2]`.
        let l_read = if p.read_nmates == 2 && loaded.len() >= 2 {
            self.read_length[0] + self.read_length[1] + 1
        } else {
            self.read_length[0]
        };
        if l_read as usize > DEF_READ_SEQ_LENGTH_MAX {
            anyhow::bail!(
                "EXITING because of FATAL ERROR in reads input: Lread = {} exceeds DEF_readSeqLengthMax = {}",
                l_read,
                DEF_READ_SEQ_LENGTH_MAX
            );
        }
        self.l_read = l_read;
        self.read_length_pair_original = if p.read_nmates == 2 && loaded.len() >= 2 {
            self.read_length_original[0] + self.read_length_original[1] + 1
        } else {
            self.read_length_original[0]
        };

        // read1[0] — numeric forward (combined)
        self.read1[0].clear();
        self.read1[0].extend_from_slice(&loaded[0].seq_num);
        if p.read_nmates == 2 && loaded.len() >= 2 {
            self.read1[0].push(MARK_FRAG_SPACER_BASE);
            let lr1 = self.read_length[1] as usize;
            let offset = self.read1[0].len();
            let mut mate2_rc = vec![0u8; lr1];
            complement_seq_numbers(&loaded[1].seq_num[..lr1], &mut mate2_rc);
            mate2_rc.reverse();
            self.read1[0].extend_from_slice(&mate2_rc);
            let _ = offset;
        } else {
            self.read_length[1] = 0;
        }

        // Ensure read1[1] / read1[2] are same length as read1[0]
        let ll = self.read1[0].len();
        self.read1[1].resize(ll, 0);
        self.read1[2].resize(ll, 0);

        // read1[1] = complement(read1[0])  (no reverse)
        let (head, tail) = self.read1.split_at_mut(1);
        complement_seq_numbers(&head[0], &mut tail[0]);
        // Re-join; subsequent lines operate on read1[2] via tail[1].
        let _ = (head, tail);

        // read1[2] = reverse(read1[1])
        for ii in 0..ll {
            self.read1[2][ll - ii - 1] = self.read1[1][ii];
        }

        // stats
        stats.read_n += 1;
        stats.read_bases += self.read_length[0] + self.read_length[1];

        // Max mismatches allowed. Port of ReadAlign_oneRead.cpp:78 —
        // C++ uses `outFilterMismatchNoverReadLmax` here (NOT `NoverLmax`).
        // In our Rust params struct that field is spelled `_nover_rmax`.
        self.out_filter_mismatch_nmax_total = std::cmp::min(
            p.out_filter_mismatch_nmax as u64,
            (p.out_filter_mismatch_nover_rmax
                * (self.read_length[0] + self.read_length[1]) as f64) as u64,
        );

        // Map
        unsafe {
            let r0 = self.read1[0].clone();
            let r1 = self.read1[1].clone();
            let r2 = self.read1[2].clone();
            let bufs: [&[u8]; 3] = [&r0, &r1, &r2];
            self.map_one_read(p, map_gen, &bufs);
        }
        // Port of `ReadAlign::peOverlapMergeMap` call site in
        // `ReadAlign_oneRead.cpp:87`. No-op unless `--peOverlapNbasesMin > 0`.
        self.pe_overlap_merge_map(p, map_gen);

        // Multimapper selection and mapped filtering. These mirror the call
        // sequence in `ReadAlign_oneRead.cpp`: multMapSelect -> mappedFilter ->
        // transcriptStats (for each selected alignment).
        self.mult_map_select(p, map_gen)?;
        let mapped = self.mapped_filter(p, stats);

        // Accumulate per-transcript stats for each selected alignment when mapped.
        if mapped {
            if self.n_tr == 1 {
                stats.mapped_reads_u += 1;
            } else if self.n_tr > 1 {
                stats.mapped_reads_m += 1;
            }
            for tr in &self.tr_mult_array {
                let exon_lengths: Vec<u64> = tr
                    .exons
                    .iter()
                    .take(tr.n_exons as usize)
                    .map(|e| e[star_core::types::EX_L])
                    .collect();
                let view = star_stats::stats::TranscriptStatView {
                    n_exons: tr.n_exons,
                    n_mm: tr.n_mm,
                    n_ins: tr.n_ins,
                    n_del: tr.n_del,
                    l_ins: tr.l_ins,
                    l_del: tr.l_del,
                    exon_lengths: &exon_lengths,
                    canon_sj: &tr.canon_sj[..tr.n_exons as usize],
                    sj_annot: &tr.sj_annot[..tr.n_exons as usize],
                };
                stats.transcript_stats(&view, self.l_read);
            }
        }

        Ok(Some(OneReadOutputs {
            mapped,
            is_chimeric: false,
        }))
    }
}
