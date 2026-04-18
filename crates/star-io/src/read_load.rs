//! 1:1 port of `readLoad.cpp` — parse one FASTA/FASTQ read from a buffered
//! stream and convert the sequence to STAR's 0..4 nucleotide encoding.

use star_core::seq::convert_nucleotides_to_numbers;
use star_params::parameters::Parameters;

/// Clip configuration carried per mate; mirrors `ClipMate` in C++ but we only
/// implement the minimum needed for M3 (no clipping yet).
#[derive(Debug, Clone, Default)]
pub struct ClipMate {
    pub clipped_info: u8,
}

impl ClipMate {
    /// `ClipMate::clip(Lread, SeqNum)` - stub; default ClipMate does nothing.
    pub fn clip(&self, _lread: u64, _seq_num: &mut [u8]) {}
}

/// Output of a successful `read_load` call.
#[derive(Debug, Clone, Default)]
pub struct LoadedRead {
    pub read_name: String,
    pub seq: Vec<u8>,
    pub seq_num: Vec<u8>,
    pub qual: Vec<u8>,
    pub lread: u64,
    pub lread_original: u64,
    pub i_read_all: u64,
    pub read_files_index: u32,
    pub read_filter: char,
    pub read_name_extra: String,
    pub file_type: i32, // 1 = FASTA, 2 = FASTQ
}

/// Port of `readLoad` (readLoad.cpp:4).
///
/// Returns `Ok(Some(record))` on a successful read, `Ok(None)` on EOF, and
/// `Err` on malformed input.
pub fn read_load(
    reader: &mut impl std::io::BufRead,
    _p: &Parameters,
    _clip_one_mate: &mut [ClipMate; 2],
) -> std::io::Result<Option<LoadedRead>> {
    let mut buf = Vec::new();

    // Line 1: header
    let n = reader.read_until(b'\n', &mut buf)?;
    if n == 0 {
        return Ok(None);
    }
    while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') {
        buf.pop();
    }
    if buf.is_empty() {
        return Ok(None);
    }
    if buf[0] != b'@' && buf[0] != b'>' {
        return Ok(None);
    }
    let header = std::str::from_utf8(&buf).map_err(io_err)?.to_string();
    let mut header_iter = header.split_whitespace();
    let read_name = header_iter.next().unwrap_or("").to_string();
    let i_read_all: u64 = header_iter.next().and_then(|s| s.parse().ok()).unwrap_or(0);
    let read_filter: char = header_iter
        .next()
        .and_then(|s| s.chars().next())
        .unwrap_or(' ');
    let read_files_index: u32 = header_iter.next().and_then(|s| s.parse().ok()).unwrap_or(0);
    let read_name_extra: String = header_iter.collect::<Vec<_>>().join(" ");

    // Line 2: sequence
    buf.clear();
    reader.read_until(b'\n', &mut buf)?;
    while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') {
        buf.pop();
    }
    let seq = buf.clone();
    let lread = seq.len() as u64;
    let mut seq_num = vec![0u8; seq.len()];
    convert_nucleotides_to_numbers(&seq, &mut seq_num);

    // FASTQ: '+' line then quality string
    let (qual, file_type) = if header.starts_with('@') {
        let mut plus = Vec::new();
        reader.read_until(b'\n', &mut plus)?;
        let mut qbuf = Vec::new();
        reader.read_until(b'\n', &mut qbuf)?;
        while qbuf.last() == Some(&b'\n') || qbuf.last() == Some(&b'\r') {
            qbuf.pop();
        }
        (qbuf, 2)
    } else {
        // FASTA: synthesize 'A' quality bytes.
        (vec![b'A'; seq.len()], 1)
    };

    Ok(Some(LoadedRead {
        read_name,
        seq,
        seq_num,
        qual,
        lread,
        lread_original: lread,
        i_read_all,
        read_files_index,
        read_filter,
        read_name_extra,
        file_type,
    }))
}

fn io_err<E: std::fmt::Display>(e: E) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fastq_record_parses() {
        let data = b"@READ1\nACGTN\n+\nIIIII\n";
        let mut rdr = std::io::BufReader::new(&data[..]);
        let p = Parameters::new();
        let mut clip = [ClipMate::default(), ClipMate::default()];
        let got = read_load(&mut rdr, &p, &mut clip).unwrap().unwrap();
        assert_eq!(got.read_name, "@READ1");
        assert_eq!(got.seq, b"ACGTN");
        assert_eq!(got.seq_num, &[0, 1, 2, 3, 4]);
        assert_eq!(got.qual, b"IIIII");
        assert_eq!(got.file_type, 2);
    }

    #[test]
    fn fasta_record_parses() {
        let data = b">r1\nACGT\n";
        let mut rdr = std::io::BufReader::new(&data[..]);
        let p = Parameters::new();
        let mut clip = [ClipMate::default(), ClipMate::default()];
        let got = read_load(&mut rdr, &p, &mut clip).unwrap().unwrap();
        assert_eq!(got.read_name, ">r1");
        assert_eq!(got.seq, b"ACGT");
        assert_eq!(got.qual, b"AAAA");
        assert_eq!(got.file_type, 1);
    }
}
