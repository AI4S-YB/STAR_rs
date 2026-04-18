//! 1:1 port of `Transcriptome.{h,cpp}` — load gene/exon info from the
//! genome directory.
//!
//! C++ source: `STAR/source/Transcriptome.cpp` (constructor, lines 7-148).
//! The constructor reads three files produced by `GTF_transcriptGeneSJ` at
//! `genomeGenerate` time (when `--sjdbGTFfile` is passed):
//! - `geneInfo.tab`
//! - `transcriptInfo.tab`   (only loaded if trSAM or GeneCounts need it)
//! - `exonInfo.tab`         (only loaded if trSAM or GeneCounts need it)
//! - `exonGeTrInfo.tab`     (only for `GeneCounts`)
//!
//! This module focuses on the GeneCounts-only subset for M7.3. The rest
//! of the structure (transcripts / exons for TranscriptomeSAM) is loaded
//! conditionally when `--quantMode TranscriptomeSAM` is active.

use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use star_params::parameters::Parameters;

/// Exon-gene structure for GeneCounts (`Transcriptome::exG` in C++).
///
/// Matches `Transcriptome.h:31-36`:
/// ```c
/// struct {
///     uint64 nEx;
///     uint64 *s, *e, *eMax;
///     uint8  *str;
///     uint32 *g, *t;
/// } exG;
/// ```
#[derive(Debug, Default, Clone)]
pub struct ExG {
    pub n_ex: u64,
    pub s: Vec<u64>,
    pub e: Vec<u64>,
    /// Running max of `e` (`e_max[i] = max(e_max[i-1], e[i])`).
    pub e_max: Vec<u64>,
    pub str_: Vec<u8>,
    pub g: Vec<u32>,
    pub t: Vec<u32>,
}

/// Port of `Transcriptome`. Currently loads the subset required for
/// `--quantMode GeneCounts`; `trS/trE/...` arrays for TranscriptomeSAM will
/// be added in M7.4.
#[derive(Debug, Default, Clone)]
pub struct Transcriptome {
    pub tr_info_dir: PathBuf,
    /// Number of genes (first line of `geneInfo.tab`).
    pub n_ge: u32,
    pub ge_id: Vec<String>,
    pub ge_name: Vec<String>,
    pub ge_biotype: Vec<String>,
    /// Exon-gene arrays for GeneCounts. Empty if GeneCounts was not requested.
    pub ex_g: ExG,
}

impl Transcriptome {
    /// Port of `Transcriptome::Transcriptome(Parameters&)` constructor
    /// (lines 7-148 of Transcriptome.cpp).
    ///
    /// Returns `Ok(None)` when `P.quant.yes == false` (no-op constructor in
    /// C++).
    pub fn load(p: &Parameters) -> Result<Option<Self>> {
        if !p.quant.yes {
            return Ok(None);
        }

        let tr_info_dir = if p.p_ge.sjdb_gtf_file == "-" {
            PathBuf::from(&p.p_ge.g_dir)
        } else {
            PathBuf::from(&p.sjdb_insert.out_dir)
        };

        // --- geneInfo.tab (always required when quant.yes) ---
        let (n_ge, ge_id, ge_name, ge_biotype) = load_gene_info(&tr_info_dir)?;

        let mut ex_g = ExG::default();
        if p.quant.ge_count.yes {
            ex_g = load_exon_ge_tr_info(&tr_info_dir)?;
        }

        Ok(Some(Self {
            tr_info_dir,
            n_ge,
            ge_id,
            ge_name,
            ge_biotype,
            ex_g,
        }))
    }
}

/// Parse `geneInfo.tab`. First line is `nGe`, next nGe lines are
/// `geID<TAB>geName<TAB>geBiotype`.
fn load_gene_info(dir: &Path) -> Result<(u32, Vec<String>, Vec<String>, Vec<String>)> {
    let path = dir.join("geneInfo.tab");
    let file = File::open(&path).with_context(|| {
        format!(
            "opening {} (GeneCounts requires --sjdbGTFfile at genomeGenerate)",
            path.display()
        )
    })?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("empty geneInfo.tab"))??;
    let n_ge: u32 = header
        .trim()
        .parse()
        .with_context(|| format!("parsing nGe from geneInfo.tab: {header:?}"))?;
    let mut ge_id = Vec::with_capacity(n_ge as usize);
    let mut ge_name = Vec::with_capacity(n_ge as usize);
    let mut ge_biotype = Vec::with_capacity(n_ge as usize);
    for _ in 0..n_ge {
        let line = lines
            .next()
            .ok_or_else(|| anyhow::anyhow!("geneInfo.tab truncated"))??;
        let mut parts = line.split_whitespace();
        let id = parts.next().unwrap_or("").to_string();
        let name = parts.next().unwrap_or("").to_string();
        let biotype = parts.next().unwrap_or("").to_string();
        ge_id.push(id);
        ge_name.push(name);
        ge_biotype.push(biotype);
    }
    Ok((n_ge, ge_id, ge_name, ge_biotype))
}

/// Parse `exonGeTrInfo.tab` (1:1 port of Transcriptome.cpp:78-98).
///
/// Line 1: `nEx`. Each subsequent line: `s e str g t` (space-separated).
/// `e_max[i]` is computed as running max of `e[0..=i]`.
fn load_exon_ge_tr_info(dir: &Path) -> Result<ExG> {
    let path = dir.join("exonGeTrInfo.tab");
    let file = File::open(&path).with_context(|| format!("opening {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    reader.read_line(&mut header)?;
    let n_ex: u64 = header
        .trim()
        .parse()
        .with_context(|| format!("parsing nEx from exonGeTrInfo.tab: {header:?}"))?;

    let mut s = Vec::with_capacity(n_ex as usize);
    let mut e = Vec::with_capacity(n_ex as usize);
    let mut e_max = Vec::with_capacity(n_ex as usize);
    let mut str_ = Vec::with_capacity(n_ex as usize);
    let mut g = Vec::with_capacity(n_ex as usize);
    let mut t = Vec::with_capacity(n_ex as usize);

    for i in 0..n_ex {
        let mut line = String::new();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            bail!("exonGeTrInfo.tab: unexpected EOF at row {i}");
        }
        let mut it = line.split_whitespace();
        let s1: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing s"))?
            .parse()?;
        let e1: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing e"))?
            .parse()?;
        let str1: i32 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing str"))?
            .parse()?;
        let g1: u32 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing g"))?
            .parse()?;
        let t1: u32 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing t"))?
            .parse()?;
        s.push(s1);
        e.push(e1);
        str_.push(str1 as u8);
        g.push(g1);
        t.push(t1);
        let prev = if i == 0 { e1 } else { *e_max.last().unwrap() };
        e_max.push(prev.max(e1));
    }

    Ok(ExG {
        n_ex,
        s,
        e,
        e_max,
        str_,
        g,
        t,
    })
}
