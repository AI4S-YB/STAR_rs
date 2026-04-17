//! Build script: populates `STAR_RS_COMPILATION_TIME_PLACE` analogous to the
//! C++ Makefile's `-D'COMPILATION_TIME_PLACE="<date> <host>:<dir>"'`.

use std::env;

fn main() {
    // Mirror the Makefile: "<iso8601-utc> <host>:<dir>" — defaulting to empty
    // strings is fine; the CLI prints this line verbatim.
    let date = env::var("SOURCE_DATE_EPOCH")
        .ok()
        .and_then(|s| s.parse::<i64>().ok())
        .map(|epoch| {
            use chrono::{DateTime, Utc};
            let dt = DateTime::<Utc>::from_timestamp(epoch, 0).unwrap_or_default();
            dt.format("%Y-%m-%dT%H:%M:%S%:z").to_string()
        })
        .unwrap_or_else(|| {
            use chrono::Utc;
            Utc::now().format("%Y-%m-%dT%H:%M:%S%:z").to_string()
        });
    let host = env::var("HOSTNAME").unwrap_or_else(|_| String::new());
    let cwd = env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| String::new());

    println!(
        "cargo:rustc-env=STAR_RS_COMPILATION_TIME_PLACE={} {}:{}",
        date, host, cwd
    );
    println!("cargo:rerun-if-env-changed=SOURCE_DATE_EPOCH");
    println!("cargo:rerun-if-env-changed=HOSTNAME");
}
