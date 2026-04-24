//! star-core: foundation layer for the STAR 1:1 Rust port.
//!
//! Ports the content of the following C++ headers / sources:
//! - `IncludeDefine.h`  -> [`types`] (constants + type aliases)
//! - `PackedArray.{h,cpp}` -> [`packed`]
//! - `SequenceFuns.{h,cpp}` -> [`seq`]
//! - `serviceFuns.cpp` (template-based) -> [`service`]
//! - `ErrorWarning.{h,cpp}` -> [`error`]
//! - `TimeFunctions.{h,cpp}` -> [`time`]
//! - `stringSubstituteAll.{h,cpp}`, `systemFunctions.{h,cpp}` -> [`util`]
//!
//! Naming discipline: `u64` for the C++ `uint` (== unsigned long long, 64-bit),
//! `i32` for `intScore`. See the `types` module for aliases.

pub mod compression;
pub mod error;
pub mod packed;
pub mod seq;
pub mod service;
pub mod time;
pub mod types;
pub mod util;

pub use types::*;
