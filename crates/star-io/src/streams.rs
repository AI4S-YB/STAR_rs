//! 1:1 port of `class InOutStreams` (InOutStreams.{h,cpp}).
//!
//! Rust representation uses buffered file handles + `Box<dyn Write>` trait
//! objects to mirror the polymorphic `ostream*` in C++. The BAM streams use
//! `noodles-bgzf`, wired once we hit M7.

use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::Path;

use parking_lot::Mutex;

/// Port of `class InOutStreams`.
pub struct InOutStreams {
    pub log_std_out: Option<Mutex<BufWriter<Box<dyn Write + Send>>>>,
    pub out_sam: Option<Mutex<BufWriter<File>>>,
    pub out_chim_sam: Option<Mutex<BufWriter<File>>>,
    pub out_chim_junction: Option<Mutex<BufWriter<File>>>,
    pub log_main: Option<Mutex<BufWriter<File>>>,
    pub log_progress: Option<Mutex<BufWriter<File>>>,
    pub log_final: Option<Mutex<BufWriter<File>>>,
    pub out_unmapped: [Option<Mutex<BufWriter<File>>>; 2],

    /// Unsorted / coord-sorted BAM streams — wired in M7. Kept as tagged
    /// option for parity with the C++ class members.
    pub out_bam_unsorted: Option<()>,
    pub out_bam_coord: Option<()>,
    pub out_quant_bam: Option<()>,
}

impl InOutStreams {
    pub fn new() -> Self {
        Self {
            log_std_out: None,
            out_sam: None,
            out_chim_sam: None,
            out_chim_junction: None,
            log_main: None,
            log_progress: None,
            log_final: None,
            out_unmapped: [None, None],
            out_bam_unsorted: None,
            out_bam_coord: None,
            out_quant_bam: None,
        }
    }

    /// Open a text file for writing (truncates existing).
    pub fn open_text(path: &Path) -> io::Result<Mutex<BufWriter<File>>> {
        let f = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)?;
        Ok(Mutex::new(BufWriter::new(f)))
    }
}

impl Default for InOutStreams {
    fn default() -> Self {
        Self::new()
    }
}
