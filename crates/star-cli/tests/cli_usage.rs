//! Integration test: `STAR --help` / `STAR` (no args) output matches the C++
//! reference byte-for-byte, modulo the compilation-time/place line.
//!
//! The reference binary path is read from `STAR_REF_BIN`. If unset, the test
//! is skipped so CI works without the upstream artefact.

use std::path::PathBuf;
use std::process::Command;

fn ref_bin() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("STAR_REF_BIN") {
        let pb = PathBuf::from(p);
        if pb.exists() {
            return Some(pb);
        }
    }
    // Fallback to the checked-in tarball binary.
    let guess = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(3)
        .map(|p| p.join("STAR/bin/Linux_x86_64/STAR"));
    guess.filter(|p| p.exists())
}

fn rust_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_star"))
}

fn filter_allowed(s: &str) -> String {
    s.lines()
        .filter(|line| !line.starts_with("STAR compilation time,server,dir="))
        .collect::<Vec<_>>()
        .join("\n")
}

fn run(bin: &PathBuf, args: &[&str]) -> String {
    let out = Command::new(bin).args(args).output().expect("exec");
    String::from_utf8_lossy(&out.stdout).into_owned()
}

#[test]
fn usage_no_args_matches_reference() {
    let Some(bin) = ref_bin() else {
        eprintln!("STAR reference binary missing, skipping");
        return;
    };
    let o = run(&bin, &[]);
    let r = run(&rust_bin(), &[]);
    assert_eq!(filter_allowed(&o), filter_allowed(&r));
}

#[test]
fn usage_help_matches_reference() {
    let Some(bin) = ref_bin() else {
        eprintln!("STAR reference binary missing, skipping");
        return;
    };
    let o = run(&bin, &["--help"]);
    let r = run(&rust_bin(), &["--help"]);
    assert_eq!(
        filter_allowed(&o),
        filter_allowed(&r),
        "STAR --help output must match byte-for-byte (modulo compile-time line)"
    );
}

#[test]
fn usage_dash_h_matches_help() {
    let Some(bin) = ref_bin() else {
        eprintln!("STAR reference binary missing, skipping");
        return;
    };
    let o = run(&bin, &["-h"]);
    let r = run(&rust_bin(), &["-h"]);
    assert_eq!(filter_allowed(&o), filter_allowed(&r));
}
