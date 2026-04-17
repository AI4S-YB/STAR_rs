//! Ports `class ReadAlignChunk` (ReadAlignChunk.{h,cpp}) for the M3
//! single-threaded SE path.
//!
//! C++ design: each worker thread owns a `ReadAlignChunk*`, which in turn owns
//! a `ReadAlign`, per-chunk output buffers, a `BAMoutput` and a
//! `ChimericDetection`. `processChunks()` is the per-thread main loop;
//! `mapChunk()` maps one block of reads.
//!
//! In Rust we make `process_chunks` a streaming function: it pulls reads
//! from one or more `BufRead` streams and writes SAM output to the provided
//! writer. The C++ "chunkIn[] / chunkOutBAM[] / istringstream" double-buffering
//! is unnecessary once we drop pthread and use `BufRead` directly.

use std::io::{BufRead, Write};

use star_genome::Genome;
use star_params::parameters::Parameters;
use star_sjdb::out_sj::OutSJ;
use star_stats::stats::Stats;

use crate::read_align::ReadAlign;

pub struct ReadAlignChunk {
    pub i_thread: u32,
    pub read_align: ReadAlign,
    pub chunk_out_sj: OutSJ,
    /// In-memory `Chimeric.out.junction` buffer (one line per chimeric read).
    /// Equivalent to `chunkOutChimJunction` in C++. Flushed by the driver.
    pub chunk_out_chim_junction: Vec<u8>,
    pub stats_ra: Stats,
    pub sam_attr_order: Vec<u16>,
}

impl ReadAlignChunk {
    pub fn new(i_thread: u32, run_rng_seed: u32, sam_attr_order: Vec<u16>) -> Self {
        Self {
            i_thread,
            read_align: ReadAlign::new(run_rng_seed.wrapping_add(i_thread)),
            chunk_out_sj: OutSJ::new(),
            chunk_out_chim_junction: Vec::new(),
            stats_ra: Stats::new(),
            sam_attr_order,
        }
    }

    /// Port of `ReadAlignChunk::processChunks` (read-map-write loop).
    ///
    /// `readers` are the per-mate input streams (size `readNends`). `sam_out`
    /// is the SAM output writer. The optional `chim_hook` is invoked once
    /// per read after `mappedFilter` and before `outputAlignments`, giving
    /// the caller a chance to run chimeric detection and append lines to
    /// `self.chunk_out_chim_junction`. Returns number of reads processed.
    pub fn process_chunks<R: BufRead, F, Q>(
        &mut self,
        p: &Parameters,
        map_gen: &Genome,
        readers: &mut [R],
        sam_out: &mut impl Write,
        mut chim_hook: F,
        mut quant_hook: Q,
    ) -> anyhow::Result<u64>
    where
        F: FnMut(&mut ReadAlign, &mut Vec<u8>) -> anyhow::Result<()>,
        Q: FnMut(&mut ReadAlign) -> anyhow::Result<()>,
    {
        self.read_align.init_from_params(p);
        self.stats_ra.reset_n();
        let mut n_reads: u64 = 0;
        loop {
            let outcome = self
                .read_align
                .one_read(p, map_gen, readers, &mut self.stats_ra)?;
            if outcome.is_none() {
                break;
            }
            self.read_align.i_read += 1;
            n_reads += 1;

            chim_hook(&mut self.read_align, &mut self.chunk_out_chim_junction)?;
            quant_hook(&mut self.read_align)?;

            self.read_align.output_alignments(
                p,
                map_gen,
                sam_out,
                &mut self.chunk_out_sj,
                &self.sam_attr_order,
            )?;
        }
        Ok(n_reads)
    }
}
