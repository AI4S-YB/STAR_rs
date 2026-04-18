//! Unit tests for the Rust STAR BAM sort algorithm (phase 2a).
//! These target algorithm-level parity with STAR's C++ (BAMoutput.cpp +
//! BAMbinSortByCoordinate.cpp + BAMbinSortUnmapped.cpp) using hand-built
//! `noodles_bam::Record`s; they do not depend on a real STAR run.

use star_bam::bam_sort_record::{read_tmp_record, write_tmp_record, BinRecordMeta, TmpRecord, UNMAPPED_SENTINEL};

#[test]
fn tmp_record_roundtrip_mapped() {
    let mut buf: Vec<u8> = Vec::new();
    let bam_bytes = vec![0xAA, 0xBB, 0xCC, 0xDD];
    let tmp = TmpRecord { ref_id: 2, pos: 1000, i_read: 7, bam_bytes: bam_bytes.clone() };
    write_tmp_record(&mut buf, &tmp).unwrap();
    assert_eq!(buf.len(), 24 + bam_bytes.len());
    let mut cursor = std::io::Cursor::new(&buf);
    let got = read_tmp_record(&mut cursor).unwrap().expect("record");
    assert_eq!(got.ref_id, 2);
    assert_eq!(got.pos, 1000);
    assert_eq!(got.i_read, 7);
    assert_eq!(got.bam_bytes, bam_bytes);
}

#[test]
fn tmp_record_roundtrip_unmapped() {
    let mut buf: Vec<u8> = Vec::new();
    let tmp = TmpRecord { ref_id: UNMAPPED_SENTINEL, pos: UNMAPPED_SENTINEL, i_read: 42, bam_bytes: vec![1, 2, 3] };
    write_tmp_record(&mut buf, &tmp).unwrap();
    let mut cursor = std::io::Cursor::new(&buf);
    let got = read_tmp_record(&mut cursor).unwrap().expect("record");
    assert_eq!(got.ref_id, UNMAPPED_SENTINEL);
    assert_eq!(got.pos, UNMAPPED_SENTINEL);
    assert_eq!(got.i_read, 42);
}

#[test]
fn tmp_record_eof_returns_none() {
    let buf: Vec<u8> = Vec::new();
    let mut cursor = std::io::Cursor::new(&buf);
    let got = read_tmp_record(&mut cursor).unwrap();
    assert!(got.is_none());
}

#[test]
fn bin_record_meta_sort_order() {
    // (coord, i_read, offset) sort: primary coord ascending, then iRead ascending.
    let mut metas = vec![
        BinRecordMeta { coord: 100, i_read: 5, file_offset: 1000 },
        BinRecordMeta { coord: 50,  i_read: 9, file_offset: 2000 },
        BinRecordMeta { coord: 100, i_read: 3, file_offset: 3000 },
        BinRecordMeta { coord: 50,  i_read: 1, file_offset: 4000 },
    ];
    metas.sort_by(|a, b| a.sort_key().cmp(&b.sort_key()));
    let coords: Vec<_> = metas.iter().map(|m| (m.coord, m.i_read)).collect();
    assert_eq!(coords, vec![(50, 1), (50, 9), (100, 3), (100, 5)]);
}

use star_bam::bam_output::BamOutput;
use std::path::PathBuf;

fn tmp_dir(test_name: &str) -> PathBuf {
    let base = std::env::temp_dir().join(format!("star-bam-test-{test_name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&base);
    std::fs::create_dir_all(&base).unwrap();
    base
}

/// Build a fake BAM record bytes payload with the noodles record's
/// serialised form. For phase-2a tests we only care about ref_id/pos
/// extraction and the bam_bytes blob round-trip.
fn fake_bam_bytes(n: usize) -> Vec<u8> {
    // 4-byte block_size prefix (value unused by our sort path) + payload.
    let mut v = vec![0u8; 4 + n];
    v[0..4].copy_from_slice(&(n as u32).to_le_bytes());
    for (i, byte) in v[4..].iter_mut().enumerate() {
        *byte = (i % 251) as u8;
    }
    v
}

#[test]
fn coord_one_align_single_bin_accumulates() {
    let dir = tmp_dir("single-bin");
    let mut out = BamOutput::new(0, &dir, /* n_bins */ 50, /* chunk_bytes */ 1 << 16).unwrap();
    // feed 3 mapped records
    for i in 0..3u64 {
        out.coord_one_align(1, 1000 + i as u32, i, fake_bam_bytes(40)).unwrap();
    }
    // all should still be in bin 0 (pre-coord_bins)
    assert_eq!(out.bin_total_n()[0], 3);
    assert_eq!(out.active_bins(), 1);
}

#[test]
fn coord_one_align_unmapped_goes_to_last_bin() {
    let dir = tmp_dir("unmapped-last-bin");
    let mut out = BamOutput::new(0, &dir, /* n_bins */ 10, /* chunk_bytes */ 1 << 16).unwrap();
    out.coord_one_align(star_bam::bam_sort_record::UNMAPPED_SENTINEL, star_bam::bam_sort_record::UNMAPPED_SENTINEL, 99, fake_bam_bytes(40)).unwrap();
    // unmapped always lands in last bin (n_bins - 1 = 9), regardless of coord_bins state
    assert_eq!(out.bin_total_n()[9], 1);
    assert_eq!(out.bin_total_n()[0], 0);
}

#[test]
fn coord_bins_splits_evenly_by_rank() {
    let dir = tmp_dir("coord-bins-even");
    // n_bins = 5 → 4 mapped bins + 1 unmapped bin.
    let mut out = BamOutput::new(0, &dir, 5, 1 << 16).unwrap();
    // Feed 8 mapped records across ref_id=1, positions 1000,2000,3000...8000.
    for i in 0..8u64 {
        out.coord_one_align(1, 1000 + 1000 * i as u32, i, fake_bam_bytes(32)).unwrap();
    }
    assert_eq!(out.active_bins(), 1);
    out.coord_bins().unwrap();
    assert_eq!(out.active_bins(), 5);
    // After coord_bins, bin 0 is the lowest-coord quarter (positions 1000,2000),
    // bin 1 is the next quarter, etc. Each of bins 0..=3 should now hold 2 records.
    let counts = out.bin_total_n();
    assert_eq!(counts[0..4].iter().sum::<u64>(), 8);
    // bin 4 is unmapped; no unmapped inputs → 0.
    assert_eq!(counts[4], 0);
}

#[test]
fn coord_bins_idempotent_after_first_call() {
    let dir = tmp_dir("coord-bins-idempotent");
    let mut out = BamOutput::new(0, &dir, 4, 1 << 16).unwrap();
    for i in 0..4u64 {
        out.coord_one_align(1, 100 * i as u32 + 100, i, fake_bam_bytes(24)).unwrap();
    }
    out.coord_bins().unwrap();
    let snapshot: Vec<u64> = out.bin_total_n().to_vec();
    out.coord_bins().unwrap();
    assert_eq!(out.bin_total_n(), snapshot.as_slice());
}
