//! 1:1 port of `TimeFunctions.{h,cpp}`.
//!
//! STAR emits timestamps in a slightly nonstandard format:
//! `strftime(..., "%b %d %H:%M:%SS", localtime(&raw))` and then erases the
//! last character. `%S` already contains seconds, so the trailing literal `S`
//! duplicates the last second digit; the erase drops it. We reproduce the
//! quirk verbatim so `Log.out` lines match.

use chrono::{DateTime, Local, TimeZone};

/// `timeMonthDayTime(rawTime)`: format and drop the last character.
pub fn time_month_day_time(raw_time: i64) -> String {
    let dt: DateTime<Local> = Local
        .timestamp_opt(raw_time, 0)
        .single()
        .unwrap_or_else(|| Local.timestamp_opt(0, 0).unwrap());
    // C++ strftime("%b %d %H:%M:%SS") then erases the last char — yielding the
    // seconds already substituted by %S. Preserve the quirk exactly.
    let mut s = dt.format("%b %d %H:%M:%SS").to_string();
    s.pop();
    s
}

/// `timeMonthDayTime()` (now).
pub fn time_month_day_time_now() -> String {
    let now = Local::now();
    let mut s = now.format("%b %d %H:%M:%SS").to_string();
    s.pop();
    s
}

/// `time(&t)` equivalent: seconds since epoch.
#[inline]
pub fn time_now() -> i64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formatter_drops_trailing_char() {
        let s = time_month_day_time_now();
        // Should look like "Apr 17 08:46:23" (no trailing S).
        assert_eq!(s.len(), 15, "unexpected format: {:?}", s);
        assert!(!s.ends_with('S'));
    }
}
