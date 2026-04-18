//! Minimum-viable SAM→BAM converter for `--outSAMtype BAM Unsorted`.
//!
//! This is NOT a byte-exact port of STAR's `ReadAlign::alignBAM`
//! (`STAR/source/ReadAlign_alignBAM.cpp`), which packs its own BAM records
//! with STAR-specific tag ordering. Instead we lean on `noodles-sam`'s
//! streaming parser + `noodles-bam`'s writer to produce a standards-
//! compliant BGZF-compressed BAM. Logical content (fields + tags) will be
//! identical to what a `samtools view -Sb` conversion of our byte-exact
//! `Aligned.out.sam` would produce.
//!
//! Byte-level equivalence with C++ STAR's BAM will come later (M7
//! byte-exact lane).

use std::fs::File;
use std::io::{BufWriter, Cursor, Read, Write};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;
use sam::alignment::io::Write as _;

/// Convert a SAM byte stream (header + body) to BGZF-compressed BAM at
/// `bam_path`. `compression` is the BGZF compression level (0..=9);
/// `None` uses the default.
pub fn convert_sam_bytes_to_bam(
    sam_bytes: &[u8],
    bam_path: &std::path::Path,
    compression: Option<u32>,
) -> Result<()> {
    let cursor = Cursor::new(sam_bytes);
    let mut reader = sam::io::Reader::new(cursor);
    let header = reader
        .read_header()
        .with_context(|| "parsing SAM header from Rust STAR output")?;

    let file =
        File::create(bam_path).with_context(|| format!("creating {}", bam_path.display()))?;
    let writer_inner: Box<dyn Write> = Box::new(BufWriter::new(file));

    let bgzf_writer = if let Some(level) = compression {
        let lvl = bgzf::io::writer::CompressionLevel::try_from(level as u8).unwrap_or_default();
        bgzf::io::writer::Builder::default()
            .set_compression_level(lvl)
            .build_from_writer(writer_inner)
    } else {
        bgzf::io::Writer::new(writer_inner)
    };
    let mut bam_writer = bam::io::Writer::from(bgzf_writer);
    bam_writer.write_header(&header)?;

    let mut record = sam::Record::default();
    loop {
        let n = reader.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        bam_writer.write_alignment_record(&header, &record)?;
    }
    bam_writer.try_finish()?;
    Ok(())
}

/// Convenience wrapper around [`convert_sam_bytes_to_bam`] that reads
/// from a file on disk.
pub fn convert_sam_file_to_bam(
    sam_path: &std::path::Path,
    bam_path: &std::path::Path,
    compression: Option<u32>,
) -> Result<()> {
    let mut buf = Vec::new();
    File::open(sam_path)
        .with_context(|| format!("opening {}", sam_path.display()))?
        .read_to_end(&mut buf)?;
    convert_sam_bytes_to_bam(&buf, bam_path, compression)
}
