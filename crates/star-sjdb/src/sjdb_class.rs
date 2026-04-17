//! 1:1 port of `SjdbClass.h` + `sjdbLoadFromStream.cpp` + `sjdbLoadFromFiles.cpp`.
//!
//! `SjdbLoci` collects splice-junction annotations from several sources
//! (saved genome, `--sjdbFileChrStartEnd`, GTF, pass-1 `SJ.out.tab`) with
//! per-entry priorities:
//!
//! | source              | priority |
//! |---------------------|----------|
//! | loaded genome sjdb  |   30     |
//! | GTF / external file |   10-20  |
//! | pass-1 novel SJ     |    0     |

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Port of `class SjdbClass`.
#[derive(Debug, Default, Clone)]
pub struct SjdbLoci {
    pub chr: Vec<String>,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    /// Strand character: `'+' | '-' | '.'`.
    pub str_: Vec<u8>,
    pub priority: Vec<u8>,
    /// Parallel to `chr` — gene ids this junction belongs to. Populated by
    /// GTF loading only.
    pub gene: Vec<BTreeSet<u64>>,
}

impl SjdbLoci {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.chr.len()
    }

    pub fn is_empty(&self) -> bool {
        self.chr.is_empty()
    }

    /// Push one record. `str_c` accepts `'1','+','2','-'`; anything else
    /// becomes `'.'`.
    fn push(&mut self, chr: String, start: u64, end: u64, str_c: u8) {
        let normalized = match str_c {
            b'1' | b'+' => b'+',
            b'2' | b'-' => b'-',
            _ => b'.',
        };
        self.chr.push(chr);
        self.start.push(start);
        self.end.push(end);
        self.str_.push(normalized);
    }

    /// Port of `sjdbLoadFromStream` (sjdbLoadFromStream.cpp:2-29).
    ///
    /// Parses whitespace-separated `chr start end str` records. Blank lines
    /// are skipped, extra columns ignored.
    pub fn load_from_stream<R: BufRead>(&mut self, reader: R) -> io::Result<()> {
        for line in reader.lines() {
            let line = line?;
            let mut it = line.split_whitespace();
            let chr = match it.next() {
                Some(c) if !c.is_empty() => c.to_string(),
                _ => continue,
            };
            let start: u64 = match it.next().and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            let end: u64 = match it.next().and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            let str_c = it.next().and_then(|s| s.bytes().next()).unwrap_or(b'.');
            self.push(chr, start, end, str_c);
        }
        Ok(())
    }

    /// Load junctions from a path, then assign `priority` to every newly
    /// added record.
    pub fn load_from_path_with_priority(
        &mut self,
        path: impl AsRef<Path>,
        priority: u8,
    ) -> io::Result<()> {
        let f = File::open(&path).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!(
                    "FATAL INPUT error, could not open input file {}: {e}",
                    path.as_ref().display()
                ),
            )
        })?;
        self.load_from_stream(BufReader::new(f))?;
        self.priority.resize(self.chr.len(), priority);
        Ok(())
    }

    /// Port of `sjdbLoadFromFiles` (sjdbLoadFromFiles.cpp:6-25). Loads every
    /// file from `--sjdbFileChrStartEnd`; a single `-` entry is the
    /// "no file" sentinel.
    pub fn load_from_files(&mut self, sjdb_files: &[String]) -> io::Result<()> {
        if sjdb_files.is_empty() || sjdb_files[0] == "-" {
            return Ok(());
        }
        for f in sjdb_files {
            self.load_from_path_with_priority(f, 10)?;
        }
        Ok(())
    }
}
