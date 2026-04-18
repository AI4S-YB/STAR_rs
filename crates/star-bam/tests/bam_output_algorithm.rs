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
