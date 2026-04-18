//! 1:1 port of `Genome::genomeLoad()` (NoSharedMemory path only).
//!
//! Ports `Genome_genomeLoad.cpp` for `pGe.gLoad=="NoSharedMemory"`. The
//! `LoadAndKeep`/`LoadAndRemove`/`LoadAndExit`/`Remove` paths are deferred
//! (M8 at the earliest; see design doc).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use star_core::packed::PackedArray;
use star_core::types::GENOME_SPACING_CHAR;
use star_params::parameters::Parameters;

use crate::genome::Genome;

/// Pre/post padding on either side of `G` when loading a genome (C++ `L=200`).
///
/// NOTE: this differs from [`crate::genome::G_OFFSET`] (100) used during
/// `genomeGenerate`. The larger padding here matches
/// `Genome_genomeLoad.cpp:244-323` which does `G1=new char[nGenome+L+L]`.
pub const LOAD_L: usize = 200;
pub const LOAD_K: u8 = 6;

/// Entry point: `Genome::genomeLoad()` (NoSharedMemory path).
pub fn genome_load(p: &mut Parameters, genome: &mut Genome) -> Result<()> {
    let g_dir = PathBuf::from(&p.p_ge.g_dir);

    // --- Read genome parameters from genomeParameters.txt ---
    let par_path = g_dir.join("genomeParameters.txt");
    let p1 = read_genome_parameters(&par_path)
        .with_context(|| format!("reading {}", par_path.display()))?;

    if p1.version_genome.is_empty() {
        bail!(
            "EXITING because of FATAL ERROR: read no value for the versionGenome \
             parameter from genomeParameters.txt file\n\
             SOLUTION: please re-generate genome from scratch with the latest version of STAR"
        );
    }
    if p1.version_genome != p.version_genome {
        bail!(
            "EXITING because of FATAL ERROR: Genome version: {} is INCOMPATIBLE \
             with running STAR version: {}",
            p1.version_genome,
            p.version_genome
        );
    }

    genome.g_strand_bit = p1.g_strand_bit;

    // --- Load chromosome info ---
    chr_info_load(&g_dir, genome)?;

    // --- Copy a few key params from on-disk genome to current Parameters ---
    p.p_ge.g_sa_index_nbases = p1.g_sa_index_nbases;
    p.p_ge.g_chr_bin_nbits = p1.g_chr_bin_nbits;
    genome.genome_chr_bin_nbases = 1u64 << p.p_ge.g_chr_bin_nbits;
    p.p_ge.g_sa_sparse_d = p1.g_sa_sparse_d;
    if !p1.g_file_sizes.is_empty() {
        p.p_ge.g_file_sizes = p1.g_file_sizes.clone();
    }

    // If sjdbOverhang was not set by user, inherit from generated genome.
    if !p.sjdb_overhang_user_set && p1.sjdb_overhang > 0 {
        p.p_ge.sjdb_overhang = p1.sjdb_overhang;
    }
    genome.sjdb_overhang = p.p_ge.sjdb_overhang;
    genome.sjdb_length = if genome.sjdb_overhang == 0 {
        0
    } else {
        genome.sjdb_overhang * 2 + 1
    };

    // --- Open Genome / SA / SAindex files, learn their sizes ---
    let (mut genome_file, n_genome) = open_stream(
        &g_dir,
        "Genome",
        p.p_ge.g_file_sizes.first().copied().unwrap_or(0),
    )?;
    let (mut sa_file, n_sa_byte) = open_stream(
        &g_dir,
        "SA",
        p.p_ge.g_file_sizes.get(1).copied().unwrap_or(0),
    )?;
    let (mut sai_file, _) = open_stream(&g_dir, "SAindex", 1)?;

    genome.n_genome = n_genome;
    genome.n_sa_byte = n_sa_byte;

    // --- Read SAindex header (gSAindexNbases + genomeSAindexStart[]) ---
    let mut u8buf = [0u8; 8];
    sai_file.read_exact(&mut u8buf).context("reading SAindex")?;
    let sai_nbases_disk = u64::from_le_bytes(u8buf);
    p.p_ge.g_sa_index_nbases = sai_nbases_disk;

    let n_entries = (p.p_ge.g_sa_index_nbases + 1) as usize;
    let mut sai_start = vec![0u64; n_entries];
    for slot in sai_start.iter_mut() {
        sai_file.read_exact(&mut u8buf).context("reading SAindex")?;
        *slot = u64::from_le_bytes(u8buf);
    }
    genome.n_sai = sai_start[p.p_ge.g_sa_index_nbases as usize];
    genome.genome_sa_index_start = sai_start;

    // --- Compute derived constants (G_strand_bit, SA bit-width, SAi bit-width) ---
    if genome.g_strand_bit == 0 {
        let mut gsb = (genome.n_genome as f64).log2().floor() as u64 + 1;
        if gsb < 32 {
            gsb = 32;
        }
        genome.g_strand_bit = gsb as u8;
    }
    genome.g_strand_mask = !(1u64 << genome.g_strand_bit);
    genome.n_sa = (genome.n_sa_byte * 8) / (genome.g_strand_bit as u64 + 1);
    genome
        .sa
        .define_bits((genome.g_strand_bit + 1) as u32, genome.n_sa);

    genome.sai_mark_n_bit = genome.g_strand_bit + 1;
    genome.sai_mark_absent_bit = genome.g_strand_bit + 2;
    genome.sai_mark_n_mask_c = 1u64 << genome.sai_mark_n_bit;
    genome.sai_mark_n_mask = !genome.sai_mark_n_mask_c;
    genome.sai_mark_absent_mask_c = 1u64 << genome.sai_mark_absent_bit;
    genome.sai_mark_absent_mask = !genome.sai_mark_absent_mask_c;
    genome
        .sai
        .define_bits((genome.g_strand_bit + 3) as u32, genome.n_sai);

    // --- Allocate buffers (G1 = nGenome + 2*L padded, SA + SAi) ---
    let n_g1 = genome.n_genome as usize + LOAD_L * 2;
    genome.n_g1_alloc = n_g1 as u64;
    genome.g = vec![GENOME_SPACING_CHAR; n_g1];

    genome.sa.allocate_array();
    genome.sai.allocate_array();

    // --- Read genome / SA / SAindex into buffers ---
    let g_start = LOAD_L;
    let g_end = g_start + genome.n_genome as usize;
    genome_file
        .read_exact(&mut genome.g[g_start..g_end])
        .context("reading Genome body")?;

    // Fill the 200-byte tails with K-1 (already done via spacing char init for
    // pre-pad; post-pad region is also spacing). We overwrite both with K-1
    // explicitly to match `G1[ii]=K-1; G[nGenome+ii]=K-1;`.
    for ii in 0..LOAD_L {
        genome.g[ii] = LOAD_K - 1;
        genome.g[g_end + ii] = LOAD_K - 1;
    }

    sa_file
        .read_exact(genome.sa.as_bytes_mut())
        .context("reading SA body")?;
    sai_file
        .read_exact(genome.sai.as_bytes_mut())
        .context("reading SAindex body")?;

    // --- chrBinFill + loadSJDB ---
    genome.chr_bin_fill();
    load_sjdb(&g_dir, genome)?;

    // --- Derive window bin parameters (Genome_genomeLoad.cpp:380-410) ---
    if p.align_intron_max == 0 && p.align_mates_gap_max == 0 {
        // No change.
    } else {
        let mate_gap = if p.align_mates_gap_max == 0 {
            1000
        } else {
            p.align_mates_gap_max
        };
        let x = std::cmp::max(std::cmp::max(4u64, p.align_intron_max), mate_gap) as f64 / 4.0;
        p.win_bin_nbits = (x.log2() + 0.5).floor() as u64;
        let floor_log2 = ((genome.n_genome / 40_000 + 1) as f64).log2() + 0.5;
        p.win_bin_nbits = std::cmp::max(p.win_bin_nbits, floor_log2.floor() as u64);
    }

    if p.win_bin_nbits > p.p_ge.g_chr_bin_nbits {
        p.win_bin_nbits = p.p_ge.g_chr_bin_nbits;
    }

    if !(p.align_intron_max == 0 && p.align_mates_gap_max == 0) {
        p.win_flank_nbins = std::cmp::max(p.align_intron_max, p.align_mates_gap_max)
            / (1u64 << p.win_bin_nbits)
            + 1;
        p.win_anchor_dist_nbins = 2 * p.win_flank_nbins;
    }

    p.win_bin_chr_nbits = p.p_ge.g_chr_bin_nbits - p.win_bin_nbits;
    p.win_bin_n = genome.n_genome / (1u64 << p.win_bin_nbits) + 1;

    Ok(())
}

/// Parsed subset of `genomeParameters.txt` — mirrors the fields STAR reads
/// back in `Genome_genomeLoad.cpp`.
#[derive(Default)]
struct GenomeParamsOnDisk {
    version_genome: String,
    g_strand_bit: u8,
    g_sa_index_nbases: u64,
    g_chr_bin_nbits: u64,
    g_sa_sparse_d: u64,
    sjdb_overhang: u64,
    g_file_sizes: Vec<u64>,
}

fn read_genome_parameters(path: &Path) -> Result<GenomeParamsOnDisk> {
    let f = File::open(path)?;
    let mut out = GenomeParamsOnDisk::default();
    for line in BufReader::new(f).lines() {
        let line = line?;
        let mut it = line.split_whitespace();
        let key = match it.next() {
            Some(k) => k,
            None => continue,
        };
        if key == "###" {
            // Comment; check for the GstrandBit marker that STAR also writes.
            if let Some(next) = it.next() {
                if next == "GstrandBit" {
                    if let Some(v) = it.next() {
                        out.g_strand_bit = v.parse().unwrap_or(0);
                    }
                }
            }
            continue;
        }
        match key {
            "versionGenome" => out.version_genome = it.next().unwrap_or("").to_string(),
            "genomeSAindexNbases" => {
                out.g_sa_index_nbases = it.next().unwrap_or("0").parse().unwrap_or(0)
            }
            "genomeChrBinNbits" => {
                out.g_chr_bin_nbits = it.next().unwrap_or("0").parse().unwrap_or(0)
            }
            "genomeSAsparseD" => out.g_sa_sparse_d = it.next().unwrap_or("0").parse().unwrap_or(0),
            "sjdbOverhang" => out.sjdb_overhang = it.next().unwrap_or("0").parse().unwrap_or(0),
            "genomeFileSizes" => {
                out.g_file_sizes = it.filter_map(|s| s.parse().ok()).collect();
            }
            _ => {}
        }
    }
    Ok(out)
}

/// Open `g_dir/name` and, if `size == 0`, seek to end to learn its length.
fn open_stream(g_dir: &Path, name: &str, size: u64) -> Result<(File, u64)> {
    let path = g_dir.join(name);
    let file = File::open(&path).with_context(|| {
        format!(
            "EXITING because of FATAL ERROR: could not open genome file: {}",
            path.display()
        )
    })?;
    let sz = if size > 0 {
        size
    } else {
        file.metadata()?.len()
    };
    Ok((file, sz))
}

/// Port of `Genome::chrInfoLoad()` (Genome.cpp:139-206).
fn chr_info_load(g_dir: &Path, genome: &mut Genome) -> Result<()> {
    // chrName.txt
    let chr_name_path = g_dir.join("chrName.txt");
    let f = File::open(&chr_name_path)
        .with_context(|| format!("could not open {}", chr_name_path.display()))?;
    for line in BufReader::new(f).lines() {
        let l = line?;
        if l.is_empty() {
            break;
        }
        genome.chr_name.push(l);
    }
    genome.n_chr_real = genome.chr_name.len() as u64;

    // chrLength.txt
    genome.chr_length = vec![0u64; genome.n_chr_real as usize];
    let f = File::open(g_dir.join("chrLength.txt"))?;
    let mut reader = BufReader::new(f);
    let mut buf = String::new();
    reader.read_to_string(&mut buf)?;
    for (slot, tok) in genome.chr_length.iter_mut().zip(buf.split_whitespace()) {
        *slot = tok.parse().unwrap_or(0);
    }

    // chrStart.txt (has nChrReal+1 entries).
    genome.chr_start = vec![0u64; (genome.n_chr_real + 1) as usize];
    let f = File::open(g_dir.join("chrStart.txt"))?;
    let mut reader = BufReader::new(f);
    let mut buf = String::new();
    reader.read_to_string(&mut buf)?;
    for (slot, tok) in genome.chr_start.iter_mut().zip(buf.split_whitespace()) {
        *slot = tok.parse().unwrap_or(0);
    }

    // Rebuild chr_name_index and keep *_all copies in sync (C++ does this
    // implicitly at struct construction; we do it here for M3 alignment use).
    genome.chr_name_index = HashMap::new();
    for (ii, nm) in genome.chr_name.iter().enumerate() {
        genome.chr_name_index.insert(nm.clone(), ii as u64);
    }
    genome.chr_name_all = genome.chr_name.clone();
    genome.chr_length_all = genome.chr_length.clone();

    Ok(())
}

/// Port of `Genome::loadSJDB` (Genome_genomeLoad.cpp:471-521).
fn load_sjdb(g_dir: &Path, genome: &mut Genome) -> Result<()> {
    // No sjdb chromosomes (nGenome == chrStart[nChrReal])
    if genome.n_genome == genome.chr_start[genome.n_chr_real as usize] {
        genome.sjdb_n = 0;
        genome.sj_g_start = genome.chr_start[genome.n_chr_real as usize] + 1;
        return Ok(());
    }

    let sjdb_path = g_dir.join("sjdbInfo.txt");
    let f = File::open(&sjdb_path)
        .with_context(|| format!("could not open {}", sjdb_path.display()))?;
    let mut reader = BufReader::new(f);
    let mut buf = String::new();
    reader.read_to_string(&mut buf)?;
    let mut tokens = buf.split_whitespace();

    let sjdb_n: u64 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
    let sjdb_overhang_on_disk: u64 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
    genome.sjdb_n = sjdb_n;
    genome.sjdb_overhang = sjdb_overhang_on_disk;

    genome.sj_chr_start = genome.n_chr_real;
    genome.sj_g_start = genome.chr_start[genome.sj_chr_start as usize];

    genome.sj_d_start = vec![0u64; sjdb_n as usize];
    genome.sj_a_start = vec![0u64; sjdb_n as usize];
    genome.sjdb_start = vec![0u64; sjdb_n as usize];
    genome.sjdb_end = vec![0u64; sjdb_n as usize];
    genome.sjdb_motif = vec![0u8; sjdb_n as usize];
    genome.sjdb_shift_left = vec![0u8; sjdb_n as usize];
    genome.sjdb_shift_right = vec![0u8; sjdb_n as usize];
    genome.sjdb_strand = vec![0u8; sjdb_n as usize];

    for ii in 0..sjdb_n as usize {
        let start: u64 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let end: u64 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let d1: u16 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let d2: u16 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let d3: u16 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let d4: u16 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);

        genome.sjdb_start[ii] = start;
        genome.sjdb_end[ii] = end;
        genome.sjdb_motif[ii] = d1 as u8;
        genome.sjdb_shift_left[ii] = d2 as u8;
        genome.sjdb_shift_right[ii] = d3 as u8;
        genome.sjdb_strand[ii] = d4 as u8;

        genome.sj_d_start[ii] = start.saturating_sub(genome.sjdb_overhang);
        genome.sj_a_start[ii] = end + 1;
        if genome.sjdb_motif[ii] == 0 {
            // non-canonical; shift back
            genome.sj_d_start[ii] += genome.sjdb_shift_left[ii] as u64;
            genome.sj_a_start[ii] += genome.sjdb_shift_left[ii] as u64;
        }
    }

    Ok(())
}

// Silence unused-import warning when the file is compiled in isolation.
#[allow(dead_code)]
fn _imports(_p: PackedArray) {}
