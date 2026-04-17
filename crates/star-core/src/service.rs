//! 1:1 port of `serviceFuns.cpp` (header-only templates).
//!
//! Provides:
//! - Comparators (as ordering functions usable with `slice::sort_by`).
//! - Generic binary-search variants (`binarySearch1/1a/1b/_leLeft/Exact`).
//! - `splitString`.
//!
//! All boundary / tie-break behavior is preserved: the original C++ goes
//! *forward* on equality for the rightmost element in `binarySearch1/1a`,
//! and *left* for `binarySearch_leLeft`.

use std::cmp::Ordering;

use crate::types::{Int32, Int64, UInt32, UInt64};

/// `sum1D`.
#[inline]
pub fn sum_1d<T: Copy + std::ops::AddAssign + Default>(a: &[T]) -> T {
    let mut s = T::default();
    for v in a {
        s += *v;
    }
    s
}

/// `funCompareNumbers<T>`
#[inline]
pub fn cmp_numbers<T: Ord>(a: &T, b: &T) -> Ordering {
    a.cmp(b)
}

/// `funCompareNumbersReverse<T>`
#[inline]
pub fn cmp_numbers_reverse<T: Ord>(a: &T, b: &T) -> Ordering {
    b.cmp(a)
}

/// `funCompareArrays<T, N>`: lexicographic compare.
#[inline]
pub fn cmp_arrays<T: Ord>(a: &[T], b: &[T]) -> Ordering {
    a.cmp(b)
}

/// `funCompareArraysShift<T, N, Shift>`: lexicographic compare starting at `shift`.
#[inline]
pub fn cmp_arrays_shift<T: Ord>(a: &[T], b: &[T], shift: usize) -> Ordering {
    a[shift..].cmp(&b[shift..])
}

/// `funCompareTypeSecondFirst<T>`: compares `a[1],b[1]` then `a[0],b[0]`.
#[inline]
pub fn cmp_type_second_first<T: Ord>(a: &[T; 2], b: &[T; 2]) -> Ordering {
    a[1].cmp(&b[1]).then_with(|| a[0].cmp(&b[0]))
}

/// `splitString(s, delim, elems)` — returns the max token length.
pub fn split_string(s: &str, delim: char) -> (usize, Vec<String>) {
    let mut elems = Vec::new();
    let mut max_l: usize = 0;
    for tok in s.split(delim) {
        max_l = max_l.max(tok.len());
        elems.push(tok.to_string());
    }
    (max_l, elems)
}

/// `binarySearch1(x, X, N) -> uint32`.
///
/// Returns index of the *rightmost* element <= x among those equal to x; or
/// `u32::MAX` if `x` is outside `[X[0], X[N-1]]`.
pub fn binary_search1<T: Ord + Copy>(x: T, xs: &[T]) -> UInt32 {
    let n = xs.len();
    if n == 0 {
        return UInt32::MAX;
    }
    if x > xs[n - 1] || x < xs[0] {
        return UInt32::MAX;
    }
    let (mut i1, mut i2) = (0usize, n - 1);
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if xs[i3] > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    while i1 < n - 1 && x == xs[i1 + 1] {
        i1 += 1;
    }
    i1 as UInt32
}

/// `binarySearch_leLeft(x, X, N, i1)`: leftmost element <= x, if any.
/// Returns `Some(i1)` on success; `None` if `x` is outside `[X[0], X[N-1]]`.
pub fn binary_search_le_left<T: Ord + Copy>(x: T, xs: &[T]) -> Option<UInt32> {
    let n = xs.len();
    if n == 0 {
        return None;
    }
    if x > xs[n - 1] || x < xs[0] {
        return None;
    }
    let (mut i1, mut i2) = (0usize, n - 1);
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if xs[i3] > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    while i1 > 0 && x == xs[i1 - 1] {
        i1 -= 1;
    }
    Some(i1 as UInt32)
}

/// `binarySearch1a(x, X, N) -> int32`: last element <= x; -1 if smaller than all.
pub fn binary_search1a<T: Ord + Copy>(x: T, xs: &[T]) -> Int32 {
    let n = xs.len();
    if n == 0 {
        return -1;
    }
    if x > xs[n - 1] {
        return (n as Int32) - 1;
    }
    if x < xs[0] {
        return -1;
    }
    let (mut i1, mut i2) = (0usize, n - 1);
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if xs[i3] > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    while i1 < n - 1 && x == xs[i1 + 1] {
        i1 += 1;
    }
    i1 as Int32
}

/// `binarySearch1b(x, X, N) -> int32`: first element >= x; -1 if > X[N-1].
/// Requires X distinct (per C++ preconditions).
pub fn binary_search1b<T: Ord + Copy>(x: T, xs: &[T]) -> Int32 {
    let n = xs.len();
    if n == 0 {
        return -1;
    }
    if x > xs[n - 1] {
        return -1;
    }
    if x <= xs[0] {
        return 0;
    }
    let (mut i1, mut i2) = (0usize, n - 1);
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if xs[i3] >= x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    i2 as Int32
}

/// `binarySearchExact(x, X, N) -> int64`; -1 if not present.
pub fn binary_search_exact<T: Ord + Copy>(x: T, xs: &[T]) -> Int64 {
    let n = xs.len();
    if n == 0 {
        return -1;
    }
    if x > xs[n - 1] || x < xs[0] {
        return -1;
    }
    let (mut i1, mut i2) = (0usize, n - 1);
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if xs[i3] >= x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    if x == xs[i2] {
        i2 as Int64
    } else if x == xs[i1] {
        i1 as Int64
    } else {
        -1
    }
}

/// 2-element uint compare: primary then secondary.
pub fn cmp_uint2(a: &[UInt64; 2], b: &[UInt64; 2]) -> Ordering {
    a.cmp(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_search_variants() {
        let xs = [1u64, 3, 3, 5, 7, 9];
        assert_eq!(binary_search1(3u64, &xs), 2); // rightmost 3
        assert_eq!(binary_search1(5u64, &xs), 3);
        assert_eq!(binary_search1(2u64, &xs), 0); // predecessor index
        assert_eq!(binary_search1(0u64, &xs), UInt32::MAX);
        assert_eq!(binary_search1(100u64, &xs), UInt32::MAX);

        assert_eq!(binary_search_le_left(3u64, &xs), Some(1)); // leftmost 3
        assert_eq!(binary_search1a(3u64, &xs), 2);
        assert_eq!(binary_search1a(0u64, &xs), -1);
        assert_eq!(binary_search1a(100u64, &xs), 5);

        let xs2 = [1u64, 3, 5, 7, 9];
        assert_eq!(binary_search1b(4u64, &xs2), 2);
        assert_eq!(binary_search1b(1u64, &xs2), 0);
        assert_eq!(binary_search1b(10u64, &xs2), -1);

        assert_eq!(binary_search_exact(5u64, &xs2), 2);
        assert_eq!(binary_search_exact(4u64, &xs2), -1);
    }

    #[test]
    fn split_string_works() {
        let (maxl, toks) = split_string("a,bb,ccc", ',');
        assert_eq!(maxl, 3);
        assert_eq!(toks, vec!["a", "bb", "ccc"]);
    }
}
