//! SAM → [`BamOutput`] adapter for phase-2a coordinate-sort pipeline.
//!
//! Reads a SAM file from disk, assigns a monotonic `i_read` counter to each
//! record, encodes it into uncompressed BAM wire format, and forwards the
//! tuple `(ref_id, pos, i_read, bam_bytes)` to
//! [`BamOutput::coord_one_align`].
//!
//! **Strategy A** is used: each SAM record is serialised to raw BAM bytes
//! (4-byte `block_size` LE prefix + record body) via a `noodles_bam`
//! writer backed by an in-memory `Vec<u8>` — no BGZF layer. Task 9
//! (`assemble_final_bam`) can reconstruct `bam::Record`s from these bytes
//! directly (e.g. `bam::io::Reader` reading from a `Cursor`).

use std::fs::File;
use std::io::{BufReader, Cursor};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write as _;

use crate::bam_output::BamOutput;
use crate::bam_sort_record::UNMAPPED_SENTINEL;

/// Feed every record in `sam_path` into `out`, returning the total count.
///
/// For each record:
/// - `i_read` is a monotonic counter starting at 0.
/// - Mapped records: `ref_id` = 0-based reference sequence index, `pos` =
///   0-based alignment start.
/// - Unmapped records: `ref_id = UNMAPPED_SENTINEL`, `pos =
///   UNMAPPED_SENTINEL`.
/// - `bam_bytes` is the BAM wire encoding: 4-byte `block_size` (LE) followed
///   by the record body. Task 9 can round-trip these via
///   `bam::io::Reader::read_record`.
pub fn feed_sam_into_bam_output(
    sam_path: &std::path::Path,
    out: &mut BamOutput,
) -> Result<u64> {
    let file = File::open(sam_path)
        .with_context(|| format!("opening SAM file {}", sam_path.display()))?;
    let mut reader = sam::io::Reader::new(BufReader::new(file));
    let header = reader
        .read_header()
        .with_context(|| format!("reading SAM header from {}", sam_path.display()))?;

    // A BAM writer backed by a plain Vec<u8> — no BGZF layer.
    // Each call to write_alignment_record appends:
    //   [4-byte block_size LE] [encoded record body]
    // to the Vec. We harvest those bytes as `bam_bytes`.
    let mut bam_writer = bam::io::Writer::from(Vec::<u8>::new());

    let mut record = sam::Record::default();
    let mut i_read: u64 = 0;

    loop {
        let n = reader
            .read_record(&mut record)
            .with_context(|| "reading SAM record")?;
        if n == 0 {
            break;
        }

        // Determine mapped / unmapped via reference_sequence_id.
        let (ref_id, pos) = match record.reference_sequence_id(&header) {
            Some(Ok(rid)) => {
                // alignment_start is 1-based; convert to 0-based.
                let pos0 = match record.alignment_start() {
                    Some(Ok(p)) => usize::from(p).saturating_sub(1) as u32,
                    _ => 0u32,
                };
                (rid as u32, pos0)
            }
            _ => (UNMAPPED_SENTINEL, UNMAPPED_SENTINEL),
        };

        // Serialise the record into raw BAM bytes.
        let raw = bam_writer.get_mut();
        raw.clear();
        bam_writer
            .write_alignment_record(&header, &record)
            .with_context(|| format!("encoding BAM record at i_read={i_read}"))?;
        let bam_bytes = bam_writer.get_ref().clone();

        out.coord_one_align(ref_id, pos, i_read, bam_bytes)
            .with_context(|| format!("coord_one_align at i_read={i_read}"))?;

        i_read += 1;
    }

    Ok(i_read)
}

/// Decode a single BAM record from the wire-format bytes produced by
/// [`feed_sam_into_bam_output`] (Strategy A).
///
/// Useful for Task 9 (`assemble_final_bam`) when reconstructing the final
/// BAM from per-bin sorted temp files.
///
/// Note: `bam::io::Reader::read_record` in noodles-bam 0.88 takes only the
/// record argument (no header).
pub fn decode_bam_bytes(bam_bytes: &[u8]) -> Result<bam::Record> {
    let mut reader = bam::io::Reader::from(Cursor::new(bam_bytes));
    let mut rec = bam::Record::default();
    reader
        .read_record(&mut rec)
        .context("decoding BAM record from bam_bytes")?;
    Ok(rec)
}
