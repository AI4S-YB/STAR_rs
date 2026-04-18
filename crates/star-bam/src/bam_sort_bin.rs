//! Per-bin mapped-record sort — phase-2a port of
//! `BAMbinSortByCoordinate.cpp`.
//!
//! For one bin index, load every per-thread tmp file for that bin into a
//! single flat buffer, build a `BinRecordMeta` array keyed by
//! `(coord, i_read, offset)`, sort, and write the records back out in
//! sorted order. Output file: `{tmp_root}/bin_{ibin:04}_sorted.rbr`.
//!
//! Output format is the same `TmpRecord` codec used everywhere in the
//! phase-2a pipeline — the final-assembly step (Task 9) decodes it and
//! writes the BAM via `noodles-bam`.

use std::fs::{remove_file, File};
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::bam_sort_record::{
    pack_coord, read_tmp_record, write_tmp_record, BinRecordMeta, TmpRecord,
};

/// Port of `BAMbinSortByCoordinate`.
///
/// * `ibin` — the bin index (0..n_bins_total-2 for mapped bins).
/// * `bin_n` — expected record count (for a sanity check only).
/// * `bin_s` — expected total bytes across all thread files (sanity check).
/// * `n_threads` — number of per-thread input files to read.
/// * `tmp_root` — the parent dir containing per-thread subdirs.
pub fn sort_mapped_bin(
    ibin: u32,
    bin_n: u64,
    bin_s: u64,
    n_threads: u32,
    tmp_root: &Path,
) -> Result<()> {
    if bin_n == 0 {
        return Ok(());
    }

    // Load every per-thread tmp file for this bin into one flat buffer.
    let mut all: Vec<u8> = Vec::with_capacity(bin_s as usize);
    for it in 0..n_threads {
        let path = tmp_root.join(it.to_string()).join(ibin.to_string());
        if !path.exists() {
            continue;
        }
        let mut f = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;
        all.extend_from_slice(&buf);
    }

    // Parse records into `(BinRecordMeta, TmpRecord)`. We keep the
    // TmpRecord around so we can emit in sorted order without a second
    // parse pass; file_offset in `BinRecordMeta` becomes the index into
    // the parallel Vec, not a byte offset (equivalent for sorting).
    let mut records: Vec<TmpRecord> = Vec::with_capacity(bin_n as usize);
    let mut metas: Vec<BinRecordMeta> = Vec::with_capacity(bin_n as usize);
    let mut cursor = std::io::Cursor::new(&all[..]);
    let mut idx: u64 = 0;
    while let Some(rec) = read_tmp_record(&mut cursor)? {
        let coord = pack_coord(rec.ref_id, rec.pos);
        metas.push(BinRecordMeta { coord, i_read: rec.i_read, file_offset: idx });
        records.push(rec);
        idx += 1;
    }
    if metas.len() as u64 != bin_n {
        bail!(
            "bin {ibin}: expected {bin_n} records, found {} on disk",
            metas.len()
        );
    }

    metas.sort_by(|a, b| a.sort_key().cmp(&b.sort_key()));

    // Emit in sorted order.
    let out_path = tmp_root.join(format!("bin_{ibin:04}_sorted.rbr"));
    let out = File::create(&out_path).with_context(|| format!("create {}", out_path.display()))?;
    let mut w = BufWriter::with_capacity(1 << 16, out);
    for m in &metas {
        write_tmp_record(&mut w, &records[m.file_offset as usize])?;
    }
    w.flush()?;

    // Clean up per-thread source files.
    for it in 0..n_threads {
        let path = tmp_root.join(it.to_string()).join(ibin.to_string());
        if path.exists() {
            let _ = remove_file(&path);
        }
    }
    Ok(())
}
