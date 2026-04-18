//! Raw FFI bindings to upstream STAR's `opal` library (SIMD Smith-Waterman).
//!
//! Only the subset actually called by STAR is exposed (`opalSearchDatabase`,
//! `opalInitSearchResult`). Callers live in `star-align::clip_cr4` (M3).

#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

#[cfg(not(opal_stub))]
mod ffi {
    use std::os::raw::{c_int, c_uchar};

    /// Matches `OpalSearchResult` in `opal.h`.
    #[repr(C)]
    pub struct OpalSearchResult {
        pub score_set: c_int,
        pub score: c_int,
        pub end_location_target: c_int,
        pub end_location_query: c_int,
        pub start_location_target: c_int,
        pub start_location_query: c_int,
        pub alignment: *mut c_uchar,
        pub alignment_length: c_int,
    }

    unsafe extern "C" {
        pub fn opalInitSearchResult(result: *mut OpalSearchResult);

        pub fn opalSearchDatabase(
            query: *const c_uchar,
            query_length: c_int,
            db: *const *const c_uchar,
            db_length: c_int,
            db_seq_lengths: *const c_int,
            gap_open: c_int,
            gap_ext: c_int,
            score_matrix: *const c_int,
            alphabet_length: c_int,
            results: *mut *mut OpalSearchResult,
            search_type: c_int,
            mode_code: c_int,
            overflow_method: c_int,
        ) -> c_int;
    }

    pub const OPAL_MODE_SW: c_int = 3;
    pub const OPAL_OVERFLOW_SIMPLE: c_int = 0;
    pub const OPAL_SEARCH_ALIGNMENT: c_int = 2;
    pub const OPAL_ERR_NO_SIMD_SUPPORT: c_int = 2;
}

#[cfg(opal_stub)]
mod ffi {
    //! Stubs when upstream opal source is not present at build time. Any call
    //! to `opalSearchDatabase` will panic; users should only hit the stub path
    //! in tooling builds that do not require ClipCR4.

    use std::os::raw::{c_int, c_uchar};

    #[repr(C)]
    pub struct OpalSearchResult {
        pub score_set: c_int,
        pub score: c_int,
        pub end_location_target: c_int,
        pub end_location_query: c_int,
        pub start_location_target: c_int,
        pub start_location_query: c_int,
        pub alignment: *mut c_uchar,
        pub alignment_length: c_int,
    }

    pub unsafe fn opalInitSearchResult(result: *mut OpalSearchResult) {
        unsafe {
            (*result).score_set = 0;
            (*result).score = 0;
            (*result).end_location_target = -1;
            (*result).end_location_query = -1;
            (*result).start_location_target = -1;
            (*result).start_location_query = -1;
            (*result).alignment = std::ptr::null_mut();
            (*result).alignment_length = 0;
        }
    }

    pub unsafe fn opalSearchDatabase(
        _query: *const c_uchar,
        _query_length: c_int,
        _db: *const *const c_uchar,
        _db_length: c_int,
        _db_seq_lengths: *const c_int,
        _gap_open: c_int,
        _gap_ext: c_int,
        _score_matrix: *const c_int,
        _alphabet_length: c_int,
        _results: *mut *mut OpalSearchResult,
        _search_type: c_int,
        _mode_code: c_int,
        _overflow_method: c_int,
    ) -> c_int {
        panic!("opal-sys was built as a stub; ClipCR4 requires opal.cpp");
    }

    pub const OPAL_MODE_SW: c_int = 3;
    pub const OPAL_OVERFLOW_SIMPLE: c_int = 0;
    pub const OPAL_SEARCH_ALIGNMENT: c_int = 2;
    pub const OPAL_ERR_NO_SIMD_SUPPORT: c_int = 2;
}

pub use ffi::*;
