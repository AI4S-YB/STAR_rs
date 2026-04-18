//! Phase-2a port of `BAMbinSortUnmapped.cpp`.
//!
//! Each per-thread unmapped bin temp file is already in ascending `iRead`
//! order (records were fed to `BamOutput::coord_one_align` in input
//! order within each thread). We perform a classic k-way merge using a
//! `BinaryHeap<Reverse<_>>` keyed on `i_read`.
//!
//! Output: `{tmp_root}/bin_{ibin:04}_sorted.rbr` — same `TmpRecord`
//! codec as the mapped bins, so the final assembly step can concat them
//! with zero format plumbing.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{remove_file, File};
use std::io::{BufReader, BufWriter};
use std::path::Path;

use anyhow::{Context, Result};

use crate::bam_sort_record::{read_tmp_record, write_tmp_record, TmpRecord};

struct HeapEntry {
    i_read: u64,
    stream_idx: usize,
    rec: TmpRecord,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool { self.i_read == other.i_read && self.stream_idx == other.stream_idx }
}
impl Eq for HeapEntry {}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> { Some(self.cmp(other)) }
}
impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Smaller i_read first, then smaller stream_idx to break ties
        // deterministically.
        self.i_read
            .cmp(&other.i_read)
            .then(self.stream_idx.cmp(&other.stream_idx))
    }
}

pub fn merge_unmapped_bin(ibin: u32, n_threads: u32, tmp_root: &Path) -> Result<()> {
    // Open each per-thread file that exists and has at least one record.
    let mut streams: Vec<BufReader<File>> = Vec::with_capacity(n_threads as usize);
    let mut paths: Vec<std::path::PathBuf> = Vec::new();
    let mut heap: BinaryHeap<Reverse<HeapEntry>> = BinaryHeap::new();

    for it in 0..n_threads {
        let path = tmp_root.join(it.to_string()).join(ibin.to_string());
        if !path.exists() {
            continue;
        }
        let f = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let mut r = BufReader::with_capacity(1 << 16, f);
        if let Some(first) = read_tmp_record(&mut r)? {
            let idx = streams.len();
            heap.push(Reverse(HeapEntry { i_read: first.i_read, stream_idx: idx, rec: first }));
            streams.push(r);
            paths.push(path);
        }
    }

    let out_path = tmp_root.join(format!("bin_{ibin:04}_sorted.rbr"));
    let out = File::create(&out_path).with_context(|| format!("create {}", out_path.display()))?;
    let mut w = BufWriter::with_capacity(1 << 16, out);

    while let Some(Reverse(entry)) = heap.pop() {
        write_tmp_record(&mut w, &entry.rec)?;
        // Pull the next record from the same stream.
        if let Some(next) = read_tmp_record(&mut streams[entry.stream_idx])? {
            heap.push(Reverse(HeapEntry { i_read: next.i_read, stream_idx: entry.stream_idx, rec: next }));
        }
    }
    drop(w);

    for p in &paths {
        let _ = remove_file(p);
    }
    Ok(())
}
