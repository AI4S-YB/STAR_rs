//! 1:1 port of `binarySearch2` variant used by `sjdbBuildIndex` and
//! `insertSeqSA` — search for a `(x, y)` pair inside two parallel sorted
//! arrays. Returns the matching index, or `-1` if not present.

/// Port of `binarySearch2` (binarySearch2.cpp) with the 2-key semantics used
/// by the SJDB pipeline: `xs` is the primary sort key, `ys` is the secondary
/// key for equal `xs` entries.
pub fn binary_search2_pair(x: u64, y: u64, xs: &[u64], ys: &[u64]) -> i64 {
    if xs.is_empty() {
        return -1;
    }
    let mut lo: i64 = 0;
    let mut hi: i64 = xs.len() as i64 - 1;
    while lo <= hi {
        let mid = (lo + hi) / 2;
        let mx = xs[mid as usize];
        if mx < x {
            lo = mid + 1;
        } else if mx > x {
            hi = mid - 1;
        } else {
            // Scan through equal xs entries to find matching y.
            let mut k = mid;
            while k >= 0 && xs[k as usize] == x {
                if ys[k as usize] == y {
                    return k;
                }
                k -= 1;
            }
            let mut k = mid + 1;
            while k < xs.len() as i64 && xs[k as usize] == x {
                if ys[k as usize] == y {
                    return k;
                }
                k += 1;
            }
            return -1;
        }
    }
    -1
}
