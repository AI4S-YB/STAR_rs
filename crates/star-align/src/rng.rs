//! MT19937 (Mersenne Twister) PRNG — required for bit-exact equivalence with
//! C++ `std::mt19937` used in `ReadAlign_multMapSelect.cpp`.
//!
//! The 32-bit MT19937 state/output sequence is well-specified; both libstdc++
//! and libc++ implement the same recurrence. Seed is `Parameters.runRNGseed`
//! (default 777).

pub struct Mt19937 {
    state: [u32; 624],
    idx: usize,
}

impl Mt19937 {
    const N: usize = 624;
    const M: usize = 397;
    const MATRIX_A: u32 = 0x9908_b0df;
    const UPPER_MASK: u32 = 0x8000_0000;
    const LOWER_MASK: u32 = 0x7fff_ffff;

    /// Seed identically to `std::mt19937 g(seed)` in libstdc++.
    pub fn new(seed: u32) -> Self {
        let mut state = [0u32; Self::N];
        state[0] = seed;
        for i in 1..Self::N {
            state[i] = (1812433253u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30)))
            .wrapping_add(i as u32);
        }
        Self {
            state,
            idx: Self::N,
        }
    }

    fn generate(&mut self) {
        for i in 0..Self::N {
            let y = (self.state[i] & Self::UPPER_MASK)
                | (self.state[(i + 1) % Self::N] & Self::LOWER_MASK);
            let mut next = self.state[(i + Self::M) % Self::N] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= Self::MATRIX_A;
            }
            self.state[i] = next;
        }
        self.idx = 0;
    }

    pub fn next_u32(&mut self) -> u32 {
        if self.idx >= Self::N {
            self.generate();
        }
        let mut y = self.state[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }
}

/// 1:1 port of `std::uniform_real_distribution<double>(0,1)` applied to an
/// MT19937 state. libstdc++ / libc++ both implement this by drawing two 32-bit
/// words and building a `[0, 1)` double; we replicate the libstdc++ recipe.
pub fn rng_uniform_real_0_to_1(rng: &mut Mt19937) -> f64 {
    // libstdc++ draws (a>>5, b>>6) and maps to [0,1).
    // Reference: bits/random.tcc `generate_canonical`.
    let a = (rng.next_u32() as u64) >> 5; // 27 bits
    let b = (rng.next_u32() as u64) >> 6; // 26 bits
    (a as f64 * 67_108_864.0 + b as f64) / 9_007_199_254_740_992.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values for `std::mt19937 g(777); g();` on libstdc++.
    /// Verified against `g++ -std=c++11` on linux x86-64.
    #[test]
    fn seed_777_matches_stdlib() {
        let mut rng = Mt19937::new(777);
        let v: Vec<u32> = (0..5).map(|_| rng.next_u32()).collect();
        assert_eq!(
            v,
            vec![655_685_735, 2_776_480_559, 1_298_611_771, 862_112_678, 266_444_375]
        );
    }
}
