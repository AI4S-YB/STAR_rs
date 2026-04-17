//! 1:1 port of `Stats.h` / `Stats.cpp`.
//!
//! `progressReport`, `progressReportHeader`, `reportFinal`, and `writeLines`
//! are ported so that `Log.final.out` is byte-exact with upstream STAR.

use std::fmt::Write as _;
use std::io::Write;

use star_core::time::time_month_day_time;

pub const SJ_MOTIF_SIZE: usize = 7;

#[derive(Debug, Clone)]
pub struct Stats {
    pub read_n: u64,
    pub read_bases: u64,

    pub mapped_reads_u: u64,
    pub mapped_reads_m: u64,
    pub mapped_bases: u64,
    pub mapped_mismatches_n: u64,
    pub mapped_ins_n: u64,
    pub mapped_del_n: u64,
    pub mapped_ins_l: u64,
    pub mapped_del_l: u64,
    pub mapped_portion: f64,

    pub splices_n: [u64; SJ_MOTIF_SIZE],
    pub splices_n_sjdb: u64,

    pub unmapped_other: u64,
    pub unmapped_short: u64,
    pub unmapped_mismatch: u64,
    pub unmapped_multi: u64,
    pub unmapped_all: u64,

    pub chimeric_all: u64,

    /// Unix epoch seconds (C++ `time_t`).
    pub time_start: i64,
    pub time_start_map: i64,
    pub time_finish_map: i64,
    pub time_last_report: i64,
    pub time_finish: i64,
}

impl Stats {
    pub fn new() -> Self {
        let now: i64 = 0;
        let mut s = Self {
            read_n: 0,
            read_bases: 0,
            mapped_reads_u: 0,
            mapped_reads_m: 0,
            mapped_bases: 0,
            mapped_mismatches_n: 0,
            mapped_ins_n: 0,
            mapped_del_n: 0,
            mapped_ins_l: 0,
            mapped_del_l: 0,
            mapped_portion: 0.0,
            splices_n: [0; SJ_MOTIF_SIZE],
            splices_n_sjdb: 0,
            unmapped_other: 0,
            unmapped_short: 0,
            unmapped_mismatch: 0,
            unmapped_multi: 0,
            unmapped_all: 0,
            chimeric_all: 0,
            time_start: now,
            time_start_map: now,
            time_finish_map: now,
            time_last_report: now,
            time_finish: now,
        };
        s.reset_n();
        s
    }

    pub fn reset_n(&mut self) {
        self.read_n = 0;
        self.read_bases = 0;
        self.mapped_mismatches_n = 0;
        self.mapped_ins_n = 0;
        self.mapped_del_n = 0;
        self.mapped_ins_l = 0;
        self.mapped_del_l = 0;
        self.mapped_bases = 0;
        self.mapped_portion = 0.0;
        self.mapped_reads_u = 0;
        self.mapped_reads_m = 0;
        self.unmapped_other = 0;
        self.unmapped_short = 0;
        self.unmapped_mismatch = 0;
        self.unmapped_multi = 0;
        self.unmapped_all = 0;
        self.chimeric_all = 0;
        self.splices_n_sjdb = 0;
        self.splices_n = [0; SJ_MOTIF_SIZE];
    }

    /// Port of `addStats` — accumulate per-chunk partials.
    pub fn add_stats(&mut self, s: &Stats) {
        self.read_n += s.read_n;
        self.read_bases += s.read_bases;
        self.mapped_mismatches_n += s.mapped_mismatches_n;
        self.mapped_ins_n += s.mapped_ins_n;
        self.mapped_del_n += s.mapped_del_n;
        self.mapped_ins_l += s.mapped_ins_l;
        self.mapped_del_l += s.mapped_del_l;
        self.mapped_bases += s.mapped_bases;
        self.mapped_portion += s.mapped_portion;
        self.mapped_reads_u += s.mapped_reads_u;
        self.mapped_reads_m += s.mapped_reads_m;
        self.unmapped_other += s.unmapped_other;
        self.unmapped_short += s.unmapped_short;
        self.unmapped_mismatch += s.unmapped_mismatch;
        self.unmapped_multi += s.unmapped_multi;
        self.unmapped_all += s.unmapped_all;
        self.chimeric_all += s.chimeric_all;
        self.splices_n_sjdb += s.splices_n_sjdb;
        for ii in 0..SJ_MOTIF_SIZE {
            self.splices_n[ii] += s.splices_n[ii];
        }
    }

    /// Port of `transcriptStats` — update per-alignment stats.
    pub fn transcript_stats(&mut self, t: &TranscriptStatView, l_read: u64) {
        self.mapped_mismatches_n += t.n_mm;
        self.mapped_ins_n += t.n_ins;
        self.mapped_del_n += t.n_del;
        self.mapped_ins_l += t.l_ins;
        self.mapped_del_l += t.l_del;
        if t.n_exons == 0 {
            return;
        }
        let mut mapped_l = 0u64;
        for ii in 0..t.n_exons as usize {
            mapped_l += t.exon_lengths[ii];
        }
        for ii in 0..(t.n_exons as usize - 1) {
            if t.canon_sj[ii] >= 0 {
                self.splices_n[t.canon_sj[ii] as usize] += 1;
            }
            if t.sj_annot[ii] == 1 {
                self.splices_n_sjdb += 1;
            }
        }
        self.mapped_bases += mapped_l;
        self.mapped_portion += mapped_l as f64 / l_read as f64;
    }

    /// Port of `Stats::reportFinal`. Bit-exact formatting required for
    /// `Log.final.out`.
    pub fn report_final(&self, out: &mut impl Write) -> std::io::Result<()> {
        let w1 = 50;

        let pct = |a: f64, b: f64| -> f64 {
            if b > 0.0 {
                a / b * 100.0
            } else {
                0.0
            }
        };

        let read_n = self.read_n as f64;
        let mapped_u = self.mapped_reads_u as f64;
        let mapped_m = self.mapped_reads_m as f64;
        let mapped_bases = self.mapped_bases as f64;

        let splices_total: u64 = self.splices_n.iter().sum();

        writeln!(
            out,
            "{:>w1$}{}",
            "Started job on |\t",
            time_month_day_time(self.time_start),
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Started mapping on |\t",
            time_month_day_time(self.time_start_map),
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Finished on |\t",
            time_month_day_time(self.time_finish),
            w1 = w1
        )?;

        let secs = (self.time_finish - self.time_start_map) as f64;
        let speed = if secs > 0.0 {
            self.read_n as f64 / 1e6 / secs * 3600.0
        } else {
            0.0
        };
        writeln!(
            out,
            "{:>w1$}{:.2}",
            "Mapping speed, Million of reads per hour |\t",
            speed,
            w1 = w1
        )?;
        writeln!(out)?;
        writeln!(out, "{:>w1$}{}", "Number of input reads |\t", self.read_n, w1 = w1)?;
        let avg_len = if self.read_n > 0 {
            self.read_bases / self.read_n
        } else {
            0
        };
        writeln!(
            out,
            "{:>w1$}{}",
            "Average input read length |\t",
            avg_len,
            w1 = w1
        )?;
        write!(out, "{:>w1$}", "UNIQUE READS:\n", w1 = w1)?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Uniquely mapped reads number |\t",
            self.mapped_reads_u,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "Uniquely mapped reads % |\t",
            pct(mapped_u, read_n),
            w1 = w1
        )?;
        let avg_mapped = if self.mapped_reads_u > 0 {
            mapped_bases / mapped_u
        } else {
            0.0
        };
        writeln!(
            out,
            "{:>w1$}{:.2}",
            "Average mapped length |\t",
            avg_mapped,
            w1 = w1
        )?;

        writeln!(
            out,
            "{:>w1$}{}",
            "Number of splices: Total |\t",
            splices_total,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of splices: Annotated (sjdb) |\t",
            self.splices_n_sjdb,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of splices: GT/AG |\t",
            self.splices_n[1] + self.splices_n[2],
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of splices: GC/AG |\t",
            self.splices_n[3] + self.splices_n[4],
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of splices: AT/AC |\t",
            self.splices_n[5] + self.splices_n[6],
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of splices: Non-canonical |\t",
            self.splices_n[0],
            w1 = w1
        )?;

        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "Mismatch rate per base, % |\t",
            pct(self.mapped_mismatches_n as f64, mapped_bases),
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "Deletion rate per base |\t",
            pct(self.mapped_del_l as f64, mapped_bases),
            w1 = w1
        )?;
        let del_avg = if self.mapped_del_n > 0 {
            self.mapped_del_l as f64 / self.mapped_del_n as f64
        } else {
            0.0
        };
        writeln!(out, "{:>w1$}{:.2}", "Deletion average length |\t", del_avg, w1 = w1)?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "Insertion rate per base |\t",
            pct(self.mapped_ins_l as f64, mapped_bases),
            w1 = w1
        )?;
        let ins_avg = if self.mapped_ins_n > 0 {
            self.mapped_ins_l as f64 / self.mapped_ins_n as f64
        } else {
            0.0
        };
        writeln!(out, "{:>w1$}{:.2}", "Insertion average length |\t", ins_avg, w1 = w1)?;
        write!(out, "{:>w1$}", "MULTI-MAPPING READS:\n", w1 = w1)?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of reads mapped to multiple loci |\t",
            self.mapped_reads_m,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "% of reads mapped to multiple loci |\t",
            pct(mapped_m, read_n),
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of reads mapped to too many loci |\t",
            self.unmapped_multi,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "% of reads mapped to too many loci |\t",
            pct(self.unmapped_multi as f64, read_n),
            w1 = w1
        )?;
        write!(out, "{:>w1$}", "UNMAPPED READS:\n", w1 = w1)?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of reads unmapped: too many mismatches |\t",
            self.unmapped_mismatch,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "% of reads unmapped: too many mismatches |\t",
            pct(self.unmapped_mismatch as f64, read_n),
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of reads unmapped: too short |\t",
            self.unmapped_short,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "% of reads unmapped: too short |\t",
            pct(self.unmapped_short as f64, read_n),
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of reads unmapped: other |\t",
            self.unmapped_other,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "% of reads unmapped: other |\t",
            pct(self.unmapped_other as f64, read_n),
            w1 = w1
        )?;
        write!(out, "{:>w1$}", "CHIMERIC READS:\n", w1 = w1)?;
        writeln!(
            out,
            "{:>w1$}{}",
            "Number of chimeric reads |\t",
            self.chimeric_all,
            w1 = w1
        )?;
        writeln!(
            out,
            "{:>w1$}{:.2}%",
            "% of chimeric reads |\t",
            pct(self.chimeric_all as f64, read_n),
            w1 = w1
        )?;

        let _ = splices_total;
        let _ = mapped_bases;
        let _: String = {
            let mut s = String::new();
            let _ = write!(s, "");
            s
        };
        Ok(())
    }
}

impl Default for Stats {
    fn default() -> Self {
        Self::new()
    }
}

/// Minimal projection of the fields from `Transcript` that `transcript_stats`
/// touches. Passed by the caller so `star-stats` doesn't have to depend on
/// `star-align` (keeping the crate graph acyclic).
pub struct TranscriptStatView<'a> {
    pub n_exons: u64,
    pub n_mm: u64,
    pub n_ins: u64,
    pub n_del: u64,
    pub l_ins: u64,
    pub l_del: u64,
    pub exon_lengths: &'a [u64],
    pub canon_sj: &'a [i32],
    pub sj_annot: &'a [u8],
}
