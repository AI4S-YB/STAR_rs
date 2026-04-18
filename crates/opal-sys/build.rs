//! Builds `opal.cpp` from the upstream STAR source tree.
//!
//! Opal is only needed for ClipCR4 (adapter clipping via SIMD Smith-Waterman)
//! and SpliceGraph (not in scope for first milestones). We build it with
//! `-mavx2` to match the upstream Makefile.

use std::path::PathBuf;

fn main() {
    println!("cargo:rustc-check-cfg=cfg(opal_stub)");
    // Upstream source lives at <workspace-root>/../STAR/source/opal.
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let opal_dir = manifest
        .parent()
        .and_then(|p| p.parent())
        .and_then(|p| p.parent())
        .map(|p| p.join("STAR/source/opal"))
        .expect("expected STAR/source/opal to exist next to the workspace");

    if !opal_dir.join("opal.cpp").exists() {
        // Skip build if the upstream tree is absent (CI on a sparse checkout).
        // Rust-side shims in `src/lib.rs` will use the `stubs` cfg.
        println!("cargo:rustc-cfg=opal_stub");
        println!("cargo:warning=opal.cpp not found, opal-sys built as stub");
        return;
    }

    let mut build = cc::Build::new();
    build
        .cpp(true)
        .file(opal_dir.join("opal.cpp"))
        .include(&opal_dir)
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-mavx2")
        .flag_if_supported("-fno-exceptions")
        .warnings(false);
    build.compile("star_opal");

    println!(
        "cargo:rerun-if-changed={}",
        opal_dir.join("opal.cpp").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        opal_dir.join("opal.h").display()
    );
}
