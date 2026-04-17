//! Raw bytes of `source/parametersDefault`, equivalent to the C++
//! `parametersDefault_xxd` produced by `xxd -i`.
//!
//! The C++ code writes the raw bytes to stdout verbatim for `STAR --help`:
//! ```cpp
//! cout.write(reinterpret_cast<char *>(parametersDefault),
//!            parametersDefault_len);
//! ```

pub static PARAMETERS_DEFAULT: &[u8] = include_bytes!("parametersDefault");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parameters_default_length() {
        // Sanity: file from upstream STAR v2.7.11b.
        assert_eq!(PARAMETERS_DEFAULT.len(), 55_490);
    }
}
