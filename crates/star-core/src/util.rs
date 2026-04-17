//! Ports small helpers: `stringSubstituteAll.cpp` and `systemFunctions.cpp`.

use std::fs;

/// `stringSubstituteAll(str, from, to)`: in-place replace-all. When `from` is
/// empty, the C++ function returns immediately (we do the same).
pub fn string_substitute_all(s: &mut String, from: &str, to: &str) {
    if from.is_empty() {
        return;
    }
    // Allocating once is equivalent; the C++ mutates in place for efficiency
    // but has the same final string contents.
    let replaced = s.replace(from, to);
    s.clear();
    s.push_str(&replaced);
}

/// `linuxProcMemory()`: read `/proc/self/status`, join `VmPeak/VmSize/VmHWM/VmRSS`
/// lines separated by `"; "`, with a trailing newline.
pub fn linux_proc_memory() -> String {
    let mut out = String::new();
    if let Ok(s) = fs::read_to_string("/proc/self/status") {
        for line in s.lines() {
            if line.starts_with("VmPeak")
                || line.starts_with("VmSize")
                || line.starts_with("VmHWM")
                || line.starts_with("VmRSS")
            {
                out.push_str(line);
                out.push_str("; ");
            }
        }
    }
    out.push('\n');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn substitute_all_preserves_tail() {
        let mut s = String::from("abcxabcxabc");
        string_substitute_all(&mut s, "x", "yx");
        assert_eq!(s, "abcyxabcyxabc");
    }

    #[test]
    fn substitute_all_noop_on_empty_from() {
        let mut s = String::from("hello");
        string_substitute_all(&mut s, "", "X");
        assert_eq!(s, "hello");
    }
}
