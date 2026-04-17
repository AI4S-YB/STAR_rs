//! Multi-threaded driver for `alignReads`.
//!
//! Ports the `runThreadN`-aware orchestration from
//! `mapThreadsSpawn.cpp` / `ReadAlignChunk_processChunks.cpp`. The C++
//! design pairs one `ReadAlignChunk` per pthread, protected by
//! `g_threadChunks.mutexInRead` + `mutexOutSAM`. Our Rust design keeps the
//! same per-thread chunk granularity but:
//!
//! - Reads are pre-loaded (serially) on the main thread into a `Vec<Mates>`
//!   where `Mates = Vec<LoadedRead>` (one entry per mate).
//! - The records are sliced into N roughly-equal chunks and handed to N
//!   worker threads through `std::thread::scope`.
//! - Each worker owns its own `ReadAlignChunk` and writes the per-chunk SAM
//!   output into an in-memory buffer.
//! - The driver concatenates the per-chunk buffers in chunk-index order and
//!   merges the `Stats` / `OutSJ` shards — this guarantees bit-exact output
//!   across thread counts when the Rust port is run against the same input.
//!
//! Matches the `--outSAMorder PairedKeepInputOrder` guarantee of the C++
//! implementation.

use std::io::{BufRead, Write};
use std::thread;

use star_genome::Genome;
use star_io::read_load::{read_load, ClipMate, LoadedRead};
use star_params::parameters::Parameters;
use star_sjdb::out_sj::OutSJ;
use star_stats::stats::Stats;

use crate::chunk::ReadAlignChunk;

/// Pre-loaded mate records for one read.
pub type Mates = Vec<LoadedRead>;

/// Port of `readChunkFastx` / the FASTX branch of
/// `ReadAlignChunk_processChunks.cpp` — single-threaded serial read of all
/// records from the input streams.
pub fn load_all_reads<R: BufRead>(
    p: &Parameters,
    readers: &mut [R],
) -> anyhow::Result<Vec<Mates>> {
    let mut out: Vec<Mates> = Vec::new();
    loop {
        let mut mates: Mates = Vec::with_capacity(readers.len());
        let mut eof = false;
        for reader in readers.iter_mut() {
            let mut clip = [ClipMate::default(), ClipMate::default()];
            match read_load(reader, p, &mut clip)? {
                Some(rec) => mates.push(rec),
                None => {
                    eof = true;
                    break;
                }
            }
        }
        if eof {
            break;
        }
        out.push(mates);
    }
    Ok(out)
}

/// Drive `alignReads` with `n_threads` workers.
///
/// `chim_hook` runs per-read inside each worker, between `mappedFilter`
/// and `outputAlignments`, appending chimeric-junction lines to a
/// per-worker buffer.
/// `quant_init` creates per-worker quantification state; `quant_hook` is
/// invoked per-read to update that state. The driver collects all
/// per-worker quant states into the returned `Vec<Q>` (one per chunk, in
/// deterministic chunk order) so the caller can merge them.
///
/// Returns `(n_reads, stats_merged, sj_merged, chim_junction_buffers, quant_states)`.
#[allow(clippy::too_many_arguments)]
pub fn run_align_multithread<F, QInit, QHook, Q>(
    p: &Parameters,
    map_gen: &Genome,
    all_reads: Vec<Mates>,
    sam_attr_order: &[u16],
    sam_out: &mut impl Write,
    run_rng_seed: u32,
    n_threads: usize,
    chim_hook: F,
    quant_init: QInit,
    quant_hook: QHook,
) -> anyhow::Result<(u64, Stats, OutSJ, Vec<Vec<u8>>, Vec<Q>)>
where
    F: Fn(&Parameters, &Genome, &mut crate::read_align::ReadAlign, &mut Vec<u8>) -> anyhow::Result<()>
        + Sync,
    Q: Send,
    QInit: Fn() -> Q + Sync,
    QHook:
        Fn(&Parameters, &Genome, &mut crate::read_align::ReadAlign, &mut Q) -> anyhow::Result<()>
            + Sync,
{
    let n_threads = n_threads.max(1);
    let total = all_reads.len();
    let chunk_size = (total + n_threads - 1) / n_threads.max(1);
    let chunks: Vec<Vec<Mates>> = all_reads
        .chunks(chunk_size.max(1))
        .map(|c| c.to_vec())
        .collect();

    let sam_attr_order_vec = sam_attr_order.to_vec();
    let chim_hook_ref = &chim_hook;
    let quant_init_ref = &quant_init;
    let quant_hook_ref = &quant_hook;

    let results: Vec<anyhow::Result<(u64, Vec<u8>, Vec<u8>, Stats, OutSJ, Q)>> =
        thread::scope(|s| {
            let mut handles = Vec::with_capacity(chunks.len());
            for (i, chunk) in chunks.into_iter().enumerate() {
                let sam_attr = sam_attr_order_vec.clone();
                let p_ref = p;
                let map_gen_ref = map_gen;
                let h: thread::ScopedJoinHandle<
                    '_,
                    anyhow::Result<(u64, Vec<u8>, Vec<u8>, Stats, OutSJ, Q)>,
                > = s.spawn(move || {
                    run_one_chunk(
                        p_ref,
                        map_gen_ref,
                        chunk,
                        &sam_attr,
                        run_rng_seed,
                        i as u32,
                        chim_hook_ref,
                        quant_init_ref,
                        quant_hook_ref,
                    )
                });
                handles.push(h);
            }
            handles
                .into_iter()
                .map(|h| h.join().unwrap_or_else(|e| Err(anyhow::anyhow!("worker panic: {e:?}"))))
                .collect()
        });

    let mut total_n: u64 = 0;
    let mut merged_stats = Stats::new();
    let mut merged_sj = OutSJ::new();
    let mut chim_bufs: Vec<Vec<u8>> = Vec::new();
    let mut quants: Vec<Q> = Vec::new();
    for r in results {
        let (n, sam_buf, chim_buf, stats, sj, q) = r?;
        total_n += n;
        merged_stats.add_stats(&stats);
        merged_sj.append(&sj);
        sam_out.write_all(&sam_buf)?;
        chim_bufs.push(chim_buf);
        quants.push(q);
    }
    Ok((total_n, merged_stats, merged_sj, chim_bufs, quants))
}

/// Run one worker: map each pre-loaded record, accumulating SAM bytes and
/// stats.
#[allow(clippy::too_many_arguments)]
fn run_one_chunk<F, QInit, QHook, Q>(
    p: &Parameters,
    map_gen: &Genome,
    reads: Vec<Mates>,
    sam_attr_order: &[u16],
    run_rng_seed: u32,
    i_thread: u32,
    chim_hook: &F,
    quant_init: &QInit,
    quant_hook: &QHook,
) -> anyhow::Result<(u64, Vec<u8>, Vec<u8>, Stats, OutSJ, Q)>
where
    F: Fn(&Parameters, &Genome, &mut crate::read_align::ReadAlign, &mut Vec<u8>) -> anyhow::Result<()>,
    QInit: Fn() -> Q,
    QHook: Fn(&Parameters, &Genome, &mut crate::read_align::ReadAlign, &mut Q) -> anyhow::Result<()>,
{
    let mut ra_chunk = ReadAlignChunk::new(
        i_thread,
        run_rng_seed.wrapping_add(i_thread),
        sam_attr_order.to_vec(),
    );
    ra_chunk.read_align.init_from_params(p);
    ra_chunk.stats_ra.reset_n();
    let mut buf: Vec<u8> = Vec::with_capacity(reads.len() * 256);
    let mut n: u64 = 0;
    let mut q: Q = quant_init();
    for mates in reads {
        let outcome =
            ra_chunk
                .read_align
                .one_read_loaded(p, map_gen, mates, &mut ra_chunk.stats_ra)?;
        if outcome.is_none() {
            break;
        }
        ra_chunk.read_align.i_read += 1;
        n += 1;
        chim_hook(
            p,
            map_gen,
            &mut ra_chunk.read_align,
            &mut ra_chunk.chunk_out_chim_junction,
        )?;
        quant_hook(p, map_gen, &mut ra_chunk.read_align, &mut q)?;
        ra_chunk.read_align.output_alignments(
            p,
            map_gen,
            &mut buf,
            &mut ra_chunk.chunk_out_sj,
            sam_attr_order,
        )?;
    }
    Ok((
        n,
        buf,
        ra_chunk.chunk_out_chim_junction,
        ra_chunk.stats_ra,
        ra_chunk.chunk_out_sj,
        q,
    ))
}
