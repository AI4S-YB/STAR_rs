//! 1:1 port of `void Genome::genomeGenerate()` (Genome_genomeGenerate.cpp:98).
//!
//! M2 status: orchestration skeleton + byte-exact pieces (FASTA scan, chr
//! info, genome parameters, raw genome file) are fully ported. The SA sort
//! driver and SAindex construction are stubbed with a clearly-marked TODO so
//! the integration test harness can exercise everything up to the SA step.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};

use star_params::parameters::Parameters;

use crate::fasta::genome_scan_fasta_files;
use crate::genome::{Genome, G_OFFSET};
use crate::io::{
    genome_parameters_write, write_chr_info, write_genome_sequence, write_sa, write_sai,
};
use crate::sa_index::genome_sa_index;
use crate::sa_sort::{sort_suffix_array, SaSortParams};

/// `void Genome::genomeGenerate()` entry point.
///
/// Returns a fresh `Genome` with state populated (chrName/Start/Length, G,
/// nGenome, nSA, GstrandBit/Mask) and writes the following files into
/// `p.p_ge.g_dir`:
///   - `Genome` (raw sequence)
///   - `chrName.txt`, `chrStart.txt`, `chrLength.txt`, `chrNameLength.txt`
///   - `genomeParameters.txt`
/// TODO (M2 tail): write `SA`, `SAindex`, `sjdbList.out.tab` once SA build is
/// ported.
pub fn genome_generate(p: &mut Parameters) -> Result<Genome> {
    let dir = Path::new(&p.p_ge.g_dir).to_path_buf();
    fs::create_dir_all(&dir)
        .with_context(|| format!("--genomeDir: could not create {}", dir.display()))?;

    let mut genome = Genome::new();
    genome.genome_chr_bin_nbases = 1u64 << p.p_ge.g_chr_bin_nbits;
    genome.sjdb_overhang = p.p_ge.sjdb_overhang;

    // Port of Genome_genomeGenerate.cpp:120-129: force sjdbOverhang=0 when no
    // sjdb source provided.
    if p.p_ge.sjdb_file_chr_start_end.first().map(String::as_str) == Some("-")
        && p.p_ge.sjdb_gtf_file == "-"
    {
        genome.sjdb_overhang = 0;
        p.p_ge.sjdb_overhang = 0;
    }

    // Pass 1: size only.
    let n_genome = genome_scan_fasta_files(
        &p.p_ge.g_fasta_files,
        genome.genome_chr_bin_nbases,
        &mut [],
        false,
        &mut genome,
        |_s| {},
    )?;

    genome.genome_sequence_allocate(n_genome);

    // Pass 2: fill sequence. Scan writes the numeric codes to the G view
    // (starts at G_OFFSET=100 inside g1). To satisfy the borrow checker we
    // temporarily swap `genome.g` out.
    let mut g_buf = std::mem::take(&mut genome.g);
    {
        let g_view = &mut g_buf[G_OFFSET..];
        genome_scan_fasta_files(
            &p.p_ge.g_fasta_files,
            genome.genome_chr_bin_nbases,
            g_view,
            true,
            &mut genome,
            |_s| {},
        )?;
    }
    genome.g = g_buf;

    let n_genome_true: u64 = genome.chr_length.iter().sum();
    if p.p_ge.g_sa_index_nbases as f64 > (n_genome_true as f64).log2() / 2.0 - 1.0 {
        eprintln!(
            "!!!!! WARNING: --genomeSAindexNbases {} is too large for the genome size={}, recommended {}",
            p.p_ge.g_sa_index_nbases,
            n_genome_true,
            ((n_genome_true as f64).log2() / 2.0 - 1.0) as i64
        );
    }

    // Write chromosome metadata files.
    write_chr_info(&dir, &genome)?;

    // Build the "-strand" half of G (G[2*N-1-i] = complement(G[i])) for SA.
    {
        let ng = genome.n_genome as usize;
        let g = &mut genome.g[G_OFFSET..];
        for ii in 0..ng {
            let src = g[ii];
            g[2 * ng - 1 - ii] = if src < 4 { 3 - src } else { src };
        }
    }

    // Count SA positions: positions where G[ii]<4 stepping by gSAsparseD.
    let sparse_d = p.p_ge.g_sa_sparse_d as usize;
    let mut n_sa: u64 = 0;
    {
        let ng = genome.n_genome as usize;
        let g = &genome.g[G_OFFSET..G_OFFSET + 2 * ng];
        for ii in (0..2 * ng).step_by(sparse_d) {
            if g[ii] < 4 {
                n_sa += 1;
            }
        }
    }
    genome.n_sa = n_sa;

    // GstrandBit = floor(log2(nGenome + limit*sjdbLength))+1, clamped to ≥32.
    let sj_len = if genome.sjdb_length == 0 {
        // At this stage sjdbLength is unknown; the C++ uses it for GstrandBit
        // reservation only. Default to 0 for genomes without sjdb so downstream
        // masks are deterministic; sjdb path (M5) will re-initialize.
        0
    } else {
        genome.sjdb_length
    };
    let span = genome.n_genome as f64 + p.limit_sjdb_insert_nsj as f64 * sj_len as f64;
    let mut g_strand_bit = (span.ln() / std::f64::consts::LN_2).floor() as i64 + 1;
    if g_strand_bit < 32 {
        g_strand_bit = 32;
    }
    genome.g_strand_bit = g_strand_bit as u8;
    genome.g_strand_mask = !(1u64 << genome.g_strand_bit);

    // Allocate SA PackedArray sized (GstrandBit+1) bits × nSA.
    genome
        .sa
        .define_bits((genome.g_strand_bit + 1) as u32, genome.n_sa);
    genome.sa.allocate_array();
    genome.n_sa_byte = genome.sa.length_byte;

    // Write genome binary (raw numeric sequence).
    write_genome_sequence(&dir, &genome)?;

    // Flip the SA-search view: `swap(G[2N-1-i], G[i])` (line 215-217). This
    // turns the buffer into "reverse-complement half | forward half" order
    // that the comparator expects.
    {
        let ng = genome.n_genome as usize;
        let g = &mut genome.g[G_OFFSET..];
        for ii in 0..ng {
            g.swap(ii, 2 * ng - 1 - ii);
        }
    }

    // Sort the suffix array.
    let sa_params = SaSortParams {
        n_genome: genome.n_genome,
        sparse_d: p.p_ge.g_sa_sparse_d,
        suffix_length_max: p.p_ge.g_suffix_length_max,
        limit_genome_generate_ram: p.limit_genome_generate_ram,
        n_g1_alloc: genome.n_g1_alloc,
        run_thread_n: p.run_thread_n.max(1) as u32,
    };
    let n_sa_built = unsafe {
        let g_ptr = genome.g.as_ptr().add(G_OFFSET);
        sort_suffix_array(g_ptr, &mut genome.sa, &sa_params, genome.g_strand_bit, &dir)?
    };
    if n_sa_built != genome.n_sa {
        anyhow::bail!(
            "SA build produced {n_sa_built} entries, expected {}",
            genome.n_sa
        );
    }

    // Un-flip G (line 343-345) so G is back in normal order before SAindex.
    {
        let ng = genome.n_genome as usize;
        let g = &mut genome.g[G_OFFSET..];
        for ii in 0..ng {
            g.swap(ii, 2 * ng - 1 - ii);
        }
    }

    // Build SAindex. Must borrow G and SA first (temporary move so we can
    // pass `&mut genome` to the builder which populates genomeSAindexStart /
    // SAiMark* fields).
    let g_buf = std::mem::take(&mut genome.g);
    let sa_pa = std::mem::take(&mut genome.sa);
    let sai = unsafe {
        let g_ptr = g_buf.as_ptr().add(G_OFFSET);
        genome_sa_index(g_ptr, &sa_pa, p.p_ge.g_sa_index_nbases as u32, &mut genome)?
    };
    genome.g = g_buf;
    genome.sa = sa_pa;

    // Write SA and SAindex files.
    write_sa(&dir, &genome.sa)?;
    write_sai(
        &dir,
        p.p_ge.g_sa_index_nbases,
        &genome.genome_sa_index_start,
        &sai,
    )?;

    // Stash SAi in genome state (used by alignReads).
    genome.sai = sai;

    // Port of Genome_genomeGenerate.cpp:368-370 — record just Genome + SA
    // byte sizes. (SAindex size is intentionally omitted; see C++.)
    p.p_ge.g_file_sizes = vec![genome.n_genome, genome.n_sa_byte];

    // Write genomeParameters.txt.
    genome_parameters_write(&dir.join("genomeParameters.txt"), p, &genome)?;

    Ok(genome)
}
