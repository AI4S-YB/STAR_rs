//! Smoke test: opal links and `opalInitSearchResult` zeros the struct.

#[cfg(not(opal_stub))]
#[test]
fn init_search_result_zeros_fields() {
    use opal_sys::OpalSearchResult;
    let mut r = OpalSearchResult {
        score_set: 42,
        score: 42,
        end_location_target: 42,
        end_location_query: 42,
        start_location_target: 42,
        start_location_query: 42,
        alignment: std::ptr::null_mut(),
        alignment_length: 42,
    };
    unsafe { opal_sys::opalInitSearchResult(&mut r as *mut _) };
    assert_eq!(r.score_set, 0);
    assert_eq!(r.end_location_target, -1);
    assert_eq!(r.alignment_length, 0);
}
