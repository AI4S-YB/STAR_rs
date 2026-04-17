//! 1:1 port of `PackedArray.{h,cpp}`.
//!
//! A bit-packed array used for `SA` / `SAinsert` / `SApass{1,2}` / `SAi`.
//! Each element occupies `word_length` bits (variable 1..=64); unaligned
//! 8-byte loads/stores access 64 bits starting from the element's byte offset
//! and mask out the relevant bits. The backing buffer is `length_byte` bytes
//! with an 8-byte slack at the end so the final element's load never reads
//! past the end.
//!
//! Storage layout matches the C++ code byte-for-byte so memory-mapped index
//! files can be pointed to directly.

use crate::types::UInt;
use std::ptr;

/// Either an owned `Vec<u8>` (C++ `allocateArray`/`deallocateArray`) or a raw
/// pointer borrowed from mmap / shared memory (C++ `pointArray`).
enum Backing {
    Owned(Vec<u8>),
    /// Non-owning pointer to externally managed memory (e.g. mmap'd file).
    /// Safety: the holder must ensure the pointed region lives at least as
    /// long as the `PackedArray`. Lifetime is erased on purpose to mirror the
    /// original C++ design. `len` is retained for bounds reasoning by callers.
    Borrowed {
        ptr: *mut u8,
        #[allow(dead_code)]
        len: usize,
    },
}

impl Clone for Backing {
    fn clone(&self) -> Self {
        match self {
            Backing::Owned(v) => Backing::Owned(v.clone()),
            // Borrowed cloning: promote to an owned copy so lifetime
            // invariants remain upheld across clones.
            Backing::Borrowed { ptr, len } => {
                let copied = unsafe { std::slice::from_raw_parts(*ptr, *len).to_vec() };
                Backing::Owned(copied)
            }
        }
    }
}

// Borrowed raw pointers don't impl Send/Sync automatically. Mirrors C++: the
// original code freely shared `PackedArray` instances across threads only for
// read access; we enforce the same invariants manually.
unsafe impl Send for Backing {}
unsafe impl Sync for Backing {}

/// 1:1 port of the C++ class.
#[derive(Clone)]
pub struct PackedArray {
    /// `wordLength`: bits per element (`N` in `defineBits`).
    pub word_length: u32,
    /// `wordCompLength = 64 - wordLength`.
    pub word_comp_length: u32,
    /// `bitRecMask = (~0u64) >> word_comp_length`.
    pub bit_rec_mask: UInt,
    /// `length`: number of elements.
    pub length: u64,
    /// `lengthByte`: size of the backing buffer in bytes.
    pub length_byte: u64,
    backing: Option<Backing>,
}

impl PackedArray {
    /// Default C++ ctor: `charArray=NULL; arrayAllocated=false`.
    pub fn new() -> Self {
        Self {
            word_length: 0,
            word_comp_length: 0,
            bit_rec_mask: 0,
            length: 0,
            length_byte: 0,
            backing: None,
        }
    }

    /// `defineBits(Nbits, lengthIn)`.
    pub fn define_bits(&mut self, n_bits: u32, length_in: u64) {
        debug_assert!(n_bits > 0 && n_bits <= 64);
        self.word_length = n_bits;
        self.word_comp_length = 64 - n_bits;
        self.bit_rec_mask = (!0u64) >> self.word_comp_length;
        self.length = length_in;
        // lengthByte = (length-1) * wordLength / 8 + sizeof(uint)
        self.length_byte = (length_in.saturating_sub(1)) * n_bits as u64 / 8 + 8;
    }

    /// `allocateArray()`: allocate `length_byte` bytes, zero the last 8.
    pub fn allocate_array(&mut self) {
        let mut buf = vec![0u8; self.length_byte as usize];
        // C++: memset(charArray+lengthByte-sizeof(uint),0,sizeof(uint));
        // Already zero since we used vec![0; ...], but keep explicit for intent.
        let tail = buf.len().saturating_sub(8);
        for b in &mut buf[tail..] {
            *b = 0;
        }
        self.backing = Some(Backing::Owned(buf));
    }

    /// `deallocateArray()`.
    pub fn deallocate_array(&mut self) {
        self.backing = None;
    }

    /// `pointArray(char*)`: borrow externally-managed memory.
    ///
    /// # Safety
    /// `ptr` must be valid for `length_byte` bytes for the lifetime of this
    /// `PackedArray`. The caller owns the memory.
    pub unsafe fn point_array(&mut self, ptr: *mut u8) {
        self.backing = Some(Backing::Borrowed {
            ptr,
            len: self.length_byte as usize,
        });
    }

    /// Raw pointer to the backing buffer (C++: `charArray`).
    #[inline]
    pub fn char_array(&self) -> *const u8 {
        match &self.backing {
            Some(Backing::Owned(v)) => v.as_ptr(),
            Some(Backing::Borrowed { ptr, .. }) => *ptr as *const u8,
            None => ptr::null(),
        }
    }

    #[inline]
    pub fn char_array_mut(&mut self) -> *mut u8 {
        match &mut self.backing {
            Some(Backing::Owned(v)) => v.as_mut_ptr(),
            Some(Backing::Borrowed { ptr, .. }) => *ptr,
            None => ptr::null_mut(),
        }
    }

    /// Byte slice view. Only safe for read-only paths that hold a `&Self`.
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        let len = self.length_byte as usize;
        if len == 0 {
            return &[];
        }
        // Safety: backing pointer is valid for `len` bytes.
        unsafe { std::slice::from_raw_parts(self.char_array(), len) }
    }

    #[inline]
    pub fn as_bytes_mut(&mut self) -> &mut [u8] {
        let len = self.length_byte as usize;
        if len == 0 {
            return &mut [];
        }
        unsafe { std::slice::from_raw_parts_mut(self.char_array_mut(), len) }
    }

    /// `PackedArray::operator[](ii)`.
    ///
    /// ```text
    /// uint b = ii * wordLength;
    /// uint B = b / 8, S = b % 8;
    /// uint a1 = *((uint*) (charArray+B));
    /// a1 = ((a1>>S) << wordCompLength) >> wordCompLength;
    /// ```
    #[inline]
    pub fn get(&self, ii: u64) -> UInt {
        let b = ii * self.word_length as u64;
        let big_b = (b / 8) as usize;
        let s = (b % 8) as u32;
        let ptr = unsafe { self.char_array().add(big_b) };
        // 8-byte unaligned load (C++ relies on x86-64 unaligned uint access).
        let a1 = unsafe { ptr::read_unaligned(ptr as *const u64) };
        ((a1 >> s) << self.word_comp_length) >> self.word_comp_length
    }

    /// `writePacked(jj, x)`.
    ///
    /// ```text
    /// b = jj*wordLength; B = b/8; S = b%8;
    /// x <<= S;
    /// *a1 = (*a1 & ~(bitRecMask<<S)) | x;
    /// ```
    #[inline]
    pub fn write_packed(&mut self, jj: u64, x: UInt) {
        let b = jj * self.word_length as u64;
        let big_b = (b / 8) as usize;
        let s = (b % 8) as u32;
        let x_shifted = x << s;
        let ptr = unsafe { self.char_array_mut().add(big_b) };
        unsafe {
            let cur = ptr::read_unaligned(ptr as *const u64);
            let new = (cur & !(self.bit_rec_mask << s)) | x_shifted;
            ptr::write_unaligned(ptr as *mut u64, new);
        }
    }
}

impl Default for PackedArray {
    fn default() -> Self {
        Self::new()
    }
}

impl std::ops::Index<u64> for PackedArray {
    type Output = u64;
    fn index(&self, _: u64) -> &u64 {
        unimplemented!("use PackedArray::get(); no long-lived reference possible for bit-packed reads")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Exhaustive round-trip for the widths STAR actually uses.
    /// A single `u64` unaligned load holds 64 bits, but `S` can be up to 7,
    /// so `wordLength + 7 <= 64` i.e. `wordLength <= 57`. STAR's SA width
    /// caps at `GstrandBit + 1 <= 40`; keep the test below that upper bound.
    #[test]
    fn roundtrip_small_widths() {
        for bits in [1u32, 2, 3, 4, 5, 7, 8, 13, 16, 29, 32, 39, 48, 55, 57] {
            let mut pa = PackedArray::new();
            let n = 64u64.min((1u64 << bits.min(16)) + 3);
            pa.define_bits(bits, n);
            pa.allocate_array();
            let mask = if bits == 64 { !0u64 } else { (1u64 << bits) - 1 };
            for i in 0..n {
                let v = (i.wrapping_mul(0x9E37_79B9_7F4A_7C15)) & mask;
                pa.write_packed(i, v);
            }
            for i in 0..n {
                let v = (i.wrapping_mul(0x9E37_79B9_7F4A_7C15)) & mask;
                assert_eq!(pa.get(i), v, "bits={} i={}", bits, i);
            }
        }
    }

    #[test]
    fn define_bits_length_byte() {
        let mut pa = PackedArray::new();
        pa.define_bits(34, 1_000_000);
        // length_byte = (999999) * 34 / 8 + 8
        assert_eq!(pa.length_byte, 999_999 * 34 / 8 + 8);
    }

    #[test]
    fn point_array_sharing() {
        let mut owner = vec![0u8; 64];
        let mut pa = PackedArray::new();
        pa.define_bits(8, 8);
        pa.length_byte = owner.len() as u64;
        unsafe { pa.point_array(owner.as_mut_ptr()) };
        pa.write_packed(3, 0x7A);
        assert_eq!(owner[3], 0x7A);
        assert_eq!(pa.get(3), 0x7A);
    }
}
