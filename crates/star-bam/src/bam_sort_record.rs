//! Temp-file record codec for the phase-2a BAM sort pipeline.
//!
//! Mirrors STAR's intent (record-payload + iRead in the same buffer) but
//! uses a Rust-native layout instead of `BAMoutput.cpp`'s raw
//! `bam_bytes || uint iRead` concatenation. Phase 2b will port the native
//! layout alongside T1's byte-exact encoder.

use std::io::{Read, Write};

use anyhow::{bail, Context, Result};

/// Sentinel value written into `ref_id`/`pos` for unmapped records
/// (matches STAR's `bamIn32[1] == (uint32)-1` check in
/// `BAMoutput::coordOneAlign`).
pub const UNMAPPED_SENTINEL: u32 = u32::MAX;

/// Magic bytes identifying the phase-2a temp-record format ("RBR1").
pub const TMP_RECORD_MAGIC: u32 = 0x5242_5231;

/// Minimum record header size on disk (magic + len + ref_id + pos + i_read).
pub const TMP_RECORD_HEADER_BYTES: usize = 24;

/// One record as it lives inside a per-thread per-bin temp file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TmpRecord {
    pub ref_id: u32,
    pub pos: u32,
    pub i_read: u64,
    pub bam_bytes: Vec<u8>,
}

/// Random-access sort key for a record already materialised inside a bin.
///
/// `file_offset` is the byte offset (from the start of the loaded bin
/// buffer) where this record's header begins, so we can rewrite records
/// in sorted order without re-parsing each one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BinRecordMeta {
    pub coord: u64,
    pub i_read: u64,
    pub file_offset: u64,
}

impl BinRecordMeta {
    /// `(coord, i_read, file_offset)` — the STAR sort key, see
    /// `BAMbinSortByCoordinate.cpp:52` (`funCompareArrays<uint,3>`).
    pub fn sort_key(&self) -> (u64, u64, u64) {
        (self.coord, self.i_read, self.file_offset)
    }
}

/// Write one `TmpRecord` into a buffered writer. Little-endian on all
/// platforms so the format is portable across dev machines.
pub fn write_tmp_record<W: Write>(w: &mut W, rec: &TmpRecord) -> Result<()> {
    let rec_len: u32 = rec
        .bam_bytes
        .len()
        .try_into()
        .context("tmp record bam_bytes exceeds u32::MAX")?;
    w.write_all(&TMP_RECORD_MAGIC.to_le_bytes())?;
    w.write_all(&rec_len.to_le_bytes())?;
    w.write_all(&rec.ref_id.to_le_bytes())?;
    w.write_all(&rec.pos.to_le_bytes())?;
    w.write_all(&rec.i_read.to_le_bytes())?;
    w.write_all(&rec.bam_bytes)?;
    Ok(())
}

/// Read one `TmpRecord`. Returns `Ok(None)` at clean EOF, `Err` on
/// truncation or magic mismatch.
pub fn read_tmp_record<R: Read>(r: &mut R) -> Result<Option<TmpRecord>> {
    let mut header = [0u8; TMP_RECORD_HEADER_BYTES];
    let mut filled = 0usize;
    while filled < TMP_RECORD_HEADER_BYTES {
        let n = r.read(&mut header[filled..])?;
        if n == 0 {
            if filled == 0 {
                return Ok(None); // clean EOF
            }
            bail!("truncated tmp record header ({filled}/24 bytes)");
        }
        filled += n;
    }
    let magic = u32::from_le_bytes(header[0..4].try_into().unwrap());
    if magic != TMP_RECORD_MAGIC {
        bail!("tmp record magic mismatch: 0x{magic:08x}");
    }
    let rec_len = u32::from_le_bytes(header[4..8].try_into().unwrap()) as usize;
    let ref_id = u32::from_le_bytes(header[8..12].try_into().unwrap());
    let pos = u32::from_le_bytes(header[12..16].try_into().unwrap());
    let i_read = u64::from_le_bytes(header[16..24].try_into().unwrap());
    let mut bam_bytes = vec![0u8; rec_len];
    r.read_exact(&mut bam_bytes).context("reading tmp record payload")?;
    Ok(Some(TmpRecord { ref_id, pos, i_read, bam_bytes }))
}

/// Helper: 64-bit coordinate = (ref_id << 32) | pos, matching STAR's
/// `alignG = (bamIn32[1]<<32) | bamIn32[2]` packing in
/// `BAMoutput::coordOneAlign`. Returns `u64::MAX` for unmapped inputs.
pub fn pack_coord(ref_id: u32, pos: u32) -> u64 {
    if ref_id == UNMAPPED_SENTINEL {
        u64::MAX
    } else {
        ((ref_id as u64) << 32) | (pos as u64)
    }
}
