//! 1:1 port of `ErrorWarning.{h,cpp}`.
//!
//! STAR's error handling logs a FATAL line to two streams, closes IO, and
//! `exit()`s with an error code. We offer two forms:
//!
//! - [`exit_with_error`] — matches the C++ semantics (stops the process).
//! - [`StarError`] — an `Err`-propagating variant used inside Rust-only call
//!   paths; the CLI layer converts it into `exit_with_error` at the top.
//!
//! Logging uses closures so the caller can route to wherever it needs
//! (`log_main`, `stderr`, etc.) without the crate pulling in higher layers.

use crate::time::time_month_day_time_now;
use std::fmt::Write as _;

#[derive(thiserror::Error, Debug)]
pub enum StarError {
    #[error("{0}")]
    Msg(String),
}

impl StarError {
    pub fn msg(s: impl Into<String>) -> Self {
        StarError::Msg(s.into())
    }
}

/// `exitWithError(message, stream1, stream2, errorInt)`.
///
/// Prints the FATAL header to the provided writers and terminates with
/// `error_int`. The C++ also deletes `P.inOut` to close files; callers in
/// Rust should drop the owning state before invoking this (or wrap in a
/// higher-level shutdown hook).
pub fn exit_with_error(
    message: &str,
    mut stream1: impl FnMut(&str),
    mut stream2: impl FnMut(&str),
    error_int: i32,
) -> ! {
    let mut buf = String::new();
    let _ = writeln!(buf, "\n{}", message);
    let _ = writeln!(
        buf,
        "{} ...... FATAL ERROR, exiting",
        time_month_day_time_now()
    );
    stream1(&buf);
    stream2(&buf);
    std::process::exit(error_int);
}

/// `warningMessage(message, stream1, stream2)`.
pub fn warning_message(
    message: &str,
    mut stream1: impl FnMut(&str),
    mut stream2: impl FnMut(&str),
) {
    let line = format!("!!!!! WARNING: {}\n", message);
    stream1(&line);
    stream2(&line);
}
