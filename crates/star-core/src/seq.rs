//! 1:1 port of `SequenceFuns.{h,cpp}`.
//!
//! Nucleotide <-> numeric conversions, reverse-complement, BAM nibble packing,
//! `chrFind`, local-alignment clippers, quality split. Semantics must be
//! byte-identical to the C++ implementation.

use crate::types::{UInt, UInt32, UInt64, MARK_FRAG_SPACER_BASE};

/// `convertNt01234(R0)`: 'A'/'a'=0, 'C'/'c'=1, 'G'/'g'=2, 'T'/'t'=3, else 4.
#[inline]
pub const fn convert_nt_01234(r0: u8) -> u8 {
    match r0 {
        b'a' | b'A' => 0,
        b'c' | b'C' => 1,
        b'g' | b'G' => 2,
        b't' | b'T' => 3,
        _ => 4,
    }
}

/// `convertNucleotidesToNumbers(R0, R1, Lread)`: in/out same length.
#[inline]
pub fn convert_nucleotides_to_numbers(r0: &[u8], r1: &mut [u8]) {
    let len = r0.len();
    debug_assert!(r1.len() >= len);
    for jj in 0..len {
        r1[jj] = convert_nt_01234(r0[jj]);
    }
}

/// In-place variant operating on a single buffer (e.g. FASTA loads).
/// Ports `convertCapitalBasesToNum` but accepts any case.
#[inline]
pub fn convert_nt_01234_inplace(buf: &mut [u8]) {
    for b in buf.iter_mut() {
        *b = convert_nt_01234(*b);
    }
}

/// `convertNucleotidesToNumbersRemoveControls(R0, R1, Lread) -> iR1`.
///
/// Copies numeric codes while **skipping control chars (< 32)**. Returns the
/// number of bytes actually written. N.B.: the C++ writes into `R1[jj]` at the
/// ORIGINAL index (not the packed `iR1`), which mirrors a latent bug we must
/// preserve for 1:1 output compatibility.
#[inline]
pub fn convert_nucleotides_to_numbers_remove_controls(r0: &[u8], r1: &mut [u8]) -> UInt {
    let mut i_r1: UInt = 0;
    for jj in 0..r0.len() {
        let c = r0[jj];
        let mapped = match c {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => {
                if (c as i32) < 32 {
                    continue;
                }
                4
            }
        };
        r1[jj] = mapped;
        i_r1 += 1;
    }
    i_r1
}

/// Complement of numeric sequences: 0<->3, 1<->2, >=4 unchanged.
#[inline]
pub fn complement_seq_numbers(reads_in: &[u8], reads_out: &mut [u8]) {
    for jj in 0..reads_in.len() {
        reads_out[jj] = match reads_in[jj] {
            3 => 0,
            2 => 1,
            1 => 2,
            0 => 3,
            x => x,
        };
    }
}

/// Reverse complement over letter-based nucleotides (A/C/G/T/N/IUPAC).
/// 1:1 port of `revComplementNucleotides(char*, char*, Lread)`.
pub fn rev_complement_nucleotides(reads_in: &[u8], reads_out: &mut [u8]) {
    let lread = reads_in.len();
    for jj in 0..lread {
        let c = reads_in[lread - 1 - jj];
        reads_out[jj] = match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'S' => b'S',
            b'W' => b'W',
            b'B' => b'V',
            b'D' => b'H',
            b'V' => b'B',
            b'H' => b'D',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'n' => b'n',
            b'r' => b'y',
            b'y' => b'r',
            b'k' => b'm',
            b'm' => b'k',
            b's' => b's',
            b'w' => b'w',
            b'b' => b'v',
            b'd' => b'h',
            b'v' => b'b',
            b'h' => b'd',
            other => other,
        };
    }
}

/// In-place std::string overload of `revComplementNucleotides`.
pub fn rev_complement_nucleotides_in_place(seq: &mut [u8]) {
    let seq1: Vec<u8> = seq.to_vec();
    let n = seq1.len();
    for jj in 0..n {
        let c = seq1[n - 1 - jj];
        seq[jj] = match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'S' => b'S',
            b'W' => b'W',
            b'B' => b'V',
            b'D' => b'H',
            b'V' => b'B',
            b'H' => b'D',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'n' => b'n',
            b'r' => b'y',
            b'y' => b'r',
            b'k' => b'm',
            b'm' => b'k',
            b's' => b's',
            b'w' => b'w',
            b'b' => b'v',
            b'd' => b'h',
            b'v' => b'b',
            b'h' => b'd',
            other => other,
        };
    }
}

/// `nuclToNumBAM`: `=ACMGRSVTWYHKDBN` -> 0..15.
#[inline]
pub const fn nucl_to_num_bam(cc: u8) -> u8 {
    match cc {
        b'=' => 0,
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'M' | b'm' => 3,
        b'G' | b'g' => 4,
        b'R' | b'r' => 5,
        b'S' | b's' => 6,
        b'V' | b'v' => 7,
        b'T' | b't' => 8,
        b'W' | b'w' => 9,
        b'Y' | b'y' => 10,
        b'H' | b'h' => 11,
        b'K' | b'k' => 12,
        b'D' | b'd' => 13,
        b'B' | b'b' => 14,
        _ => 15,
    }
}

/// Pack two nucleotides per byte, high-nibble first; odd tail gets low nibble 0.
pub fn nucl_pack_bam(reads_in: &[u8], reads_out: &mut [u8]) {
    let lread = reads_in.len();
    for jj in 0..(lread / 2) {
        reads_out[jj] =
            (nucl_to_num_bam(reads_in[2 * jj]) << 4) | nucl_to_num_bam(reads_in[2 * jj + 1]);
    }
    if lread % 2 == 1 {
        reads_out[lread / 2] = nucl_to_num_bam(reads_in[lread - 1]) << 4;
    }
}

/// `convertNuclStrToInt32`: 2-bit encode of string of length <=16.
/// Returns Ok((int_out, pos_n)) with `pos_n = -1` if no N, `pos_n = index` if
/// single N, and Err(()) if two or more Ns (C++ returns -2).
pub fn convert_nucl_str_to_int32(s: &[u8]) -> Result<(UInt32, i32), ()> {
    let mut int_out: UInt32 = 0;
    let mut pos_n: i32 = -1;
    for (ii, &c) in s.iter().enumerate() {
        let mut nt = convert_nt_01234(c) as UInt32;
        if nt > 3 {
            if pos_n >= 0 {
                return Err(());
            }
            pos_n = ii as i32;
            nt = 0;
        }
        int_out <<= 2;
        int_out += nt;
    }
    Ok((int_out, pos_n))
}

pub fn convert_nucl_int32_to_string(mut nucl_num: UInt32, l: UInt32) -> String {
    let mut out = vec![b'N'; l as usize];
    for ii in 1..=l {
        out[(l - ii) as usize] = match nucl_num & 3 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        nucl_num >>= 2;
    }
    // Safe: all bytes are ASCII A/C/G/T/N
    unsafe { String::from_utf8_unchecked(out) }
}

pub fn convert_nucl_str_to_int64(s: &[u8]) -> Result<(UInt64, i64), ()> {
    let mut int_out: UInt64 = 0;
    let mut pos_n: i64 = -1;
    for (ii, &c) in s.iter().enumerate() {
        let mut nt = convert_nt_01234(c) as UInt64;
        if nt > 3 {
            if pos_n >= 0 {
                return Err(());
            }
            pos_n = ii as i64;
            nt = 0;
        }
        int_out <<= 2;
        int_out += nt;
    }
    Ok((int_out, pos_n))
}

pub fn convert_nucl_int64_to_string(mut nucl_num: UInt64, l: UInt32) -> String {
    let mut out = vec![b'N'; l as usize];
    for ii in 1..=l {
        out[(l - ii) as usize] = match nucl_num & 3 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        nucl_num >>= 2;
    }
    unsafe { String::from_utf8_unchecked(out) }
}

/// `chrFind(Start, i2, chrStart)`: binary-search predecessor in `chrStart[0..i2]`.
pub fn chr_find(start: UInt, mut i2: UInt, chr_start: &[UInt]) -> UInt {
    let mut i1: UInt = 0;
    while i1 + 1 < i2 {
        let i3 = (i1 + i2) / 2;
        if chr_start[i3 as usize] > start {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    i1
}

/// `localSearch`: find best alignment offset in x (length nx) for y (length ny).
/// Ns in x skipped (`x[ix+iy]>3 continue`), mismatches counted.
pub fn local_search(x: &[u8], y: &[u8], p_mm: f64) -> UInt {
    let nx = x.len() as UInt;
    let ny = y.len() as UInt;
    let mut n_match_best: UInt = 0;
    let mut n_mm_best: UInt = 0;
    let mut ix_best: UInt = nx;
    for ix in 0..nx {
        let mut n_match: UInt = 0;
        let mut n_mm: UInt = 0;
        let lim = ny.min(nx - ix);
        for iy in 0..lim {
            let xv = x[(ix + iy) as usize];
            if xv > 3 {
                continue;
            }
            if xv == y[iy as usize] {
                n_match += 1;
            } else {
                n_mm += 1;
            }
        }
        let cond_better = n_match > n_match_best || (n_match == n_match_best && n_mm < n_mm_best);
        if cond_better && n_match > 0 && (n_mm as f64) / (n_match as f64) <= p_mm {
            ix_best = ix;
            n_match_best = n_match;
            n_mm_best = n_mm;
        }
    }
    ix_best
}

/// `localSearchNisMM`: like `local_search` but Ns in x or y count as MM.
pub fn local_search_n_is_mm(x: &[u8], y: &[u8], p_mm: f64) -> UInt {
    let nx = x.len() as UInt;
    let ny = y.len() as UInt;
    let mut n_match_best: UInt = 0;
    let mut n_mm_best: UInt = 0;
    let mut ix_best: UInt = nx;
    for ix in 0..nx {
        let mut n_match: UInt = 0;
        let mut n_mm: UInt = 0;
        let lim = ny.min(nx - ix);
        for iy in 0..lim {
            if x[(ix + iy) as usize] == y[iy as usize] && y[iy as usize] < 4 {
                n_match += 1;
            } else {
                n_mm += 1;
            }
        }
        let cond_better = n_match > n_match_best || (n_match == n_match_best && n_mm < n_mm_best);
        if cond_better && n_match > 0 && (n_mm as f64) / (n_match as f64) <= p_mm {
            ix_best = ix;
            n_match_best = n_match;
            n_mm_best = n_mm;
        }
    }
    ix_best
}

/// `localAlignHammingDist(text, query, pos) -> distBest` (via return value +
/// mutable `pos`). Treats 'N' in query as a free pass (not a mismatch).
pub fn local_align_hamming_dist(text: &[u8], query: &[u8]) -> (UInt32, UInt32) {
    let mut pos: UInt32 = 0;
    let q_size = query.len();
    let t_size = text.len();
    if t_size < q_size {
        return ((t_size + 1) as UInt32, pos);
    }
    let mut dist_best: UInt32 = q_size as UInt32;
    for ii in 0..=(t_size - q_size) {
        let mut dist1: UInt32 = 0;
        for jj in 0..q_size {
            if query[jj] != b'N' && text[jj + ii] != query[jj] {
                dist1 += 1;
            }
        }
        if dist1 < dist_best {
            dist_best = dist1;
            pos = ii as UInt32;
        }
    }
    (dist_best, pos)
}

/// `qualitySplit`: split a numeric read by bad bases (>3); outputs up to
/// `max_n_split` regions via `split_r[0/1/2][iS]`. Returns number of regions.
/// The three rows are `{start, length, frag}` in the C++ 2D array.
pub fn quality_split(
    r: &[u8],
    max_n_split: UInt,
    min_l_split: UInt,
    split_r: &mut [[UInt; 3]],
) -> UInt {
    let l = r.len() as UInt;
    let mut i_r: UInt = 0;
    let mut i_s: UInt = 0;
    let mut l_good_min: UInt = 0;
    let mut i_frag: UInt = 0;
    while i_r < l && i_s < max_n_split {
        while i_r < l && r[i_r as usize] > 3 {
            if r[i_r as usize] == MARK_FRAG_SPACER_BASE {
                i_frag += 1;
            }
            i_r += 1;
        }
        if i_r == l {
            break;
        }
        let i_r1 = i_r;
        while i_r < l && r[i_r as usize] <= 3 {
            i_r += 1;
        }
        if (i_r - i_r1) > l_good_min {
            l_good_min = i_r - i_r1;
        }
        if (i_r - i_r1) < min_l_split {
            continue;
        }
        split_r[i_s as usize][0] = i_r1;
        split_r[i_s as usize][1] = i_r - i_r1;
        split_r[i_s as usize][2] = i_frag;
        i_s += 1;
    }
    if i_s == 0 && !split_r.is_empty() {
        split_r[0][1] = l_good_min;
    }
    i_s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nt_01234_basic() {
        assert_eq!(convert_nt_01234(b'A'), 0);
        assert_eq!(convert_nt_01234(b'c'), 1);
        assert_eq!(convert_nt_01234(b'G'), 2);
        assert_eq!(convert_nt_01234(b't'), 3);
        assert_eq!(convert_nt_01234(b'N'), 4);
        assert_eq!(convert_nt_01234(b'X'), 4);
    }

    #[test]
    fn rev_complement_letters() {
        let input = b"ACGTN";
        let mut out = vec![0u8; input.len()];
        rev_complement_nucleotides(input, &mut out);
        assert_eq!(&out, b"NACGT");
    }

    #[test]
    fn nucl_pack_bam_works() {
        let input = b"ACGT";
        let mut out = vec![0u8; 2];
        nucl_pack_bam(input, &mut out);
        assert_eq!(out, vec![0x12, 0x48]); // A=1,C=2,G=4,T=8
        let odd = b"A";
        let mut out2 = vec![0u8; 1];
        nucl_pack_bam(odd, &mut out2);
        assert_eq!(out2, vec![0x10]);
    }

    #[test]
    fn chr_find_matches_binary_search() {
        let starts = vec![0u64, 1000, 3000, 7000];
        assert_eq!(chr_find(500, starts.len() as u64, &starts), 0);
        assert_eq!(chr_find(1500, starts.len() as u64, &starts), 1);
        assert_eq!(chr_find(3000, starts.len() as u64, &starts), 2);
        assert_eq!(chr_find(9999, starts.len() as u64, &starts), 3);
    }

    #[test]
    fn nucl_int32_roundtrip() {
        let (n, pos) = convert_nucl_str_to_int32(b"ACGT").unwrap();
        assert_eq!(pos, -1);
        assert_eq!(convert_nucl_int32_to_string(n, 4), "ACGT");
    }
}
