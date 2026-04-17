//! Chromosome info writers and `genomeParametersWrite` port.
//!
//! Ports:
//! - `Genome::writeChrInfo` (Genome_genomeGenerate.cpp:417-432) — writes
//!   `chrName.txt`, `chrStart.txt`, `chrLength.txt`, `chrNameLength.txt`.
//! - `Genome::writeGenomeSequence` — writes the raw `Genome` binary.
//! - `genomeParametersWrite.cpp` — writes `genomeParameters.txt`.
//!
//! All three emit byte-for-byte the same as upstream for the same inputs.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use star_core::packed::PackedArray;
use star_params::parameters::Parameters;

use crate::genome::Genome;

/// `Genome::writeChrInfo(dirOut)`. Writes all four chromosome metadata files.
pub fn write_chr_info(dir_out: &Path, genome: &Genome) -> Result<()> {
    let mut chr_n = open_w(&dir_out.join("chrName.txt"))?;
    let mut chr_s = open_w(&dir_out.join("chrStart.txt"))?;
    let mut chr_l = open_w(&dir_out.join("chrLength.txt"))?;
    let mut chr_nl = open_w(&dir_out.join("chrNameLength.txt"))?;

    for ii in 0..genome.n_chr_real as usize {
        writeln!(chr_n, "{}", genome.chr_name[ii])?;
        writeln!(chr_s, "{}", genome.chr_start[ii])?;
        writeln!(chr_l, "{}", genome.chr_length[ii])?;
        writeln!(chr_nl, "{}\t{}", genome.chr_name[ii], genome.chr_length[ii])?;
    }
    // Trailing end-of-genome marker in chrStart.txt.
    writeln!(chr_s, "{}", genome.chr_start[genome.n_chr_real as usize])?;
    Ok(())
}

/// `Genome::writeGenomeSequence(dirOut)`. Writes `Genome` (raw numeric bytes).
pub fn write_genome_sequence(dir_out: &Path, genome: &Genome) -> Result<()> {
    let mut f = open_w(&dir_out.join("Genome"))?;
    let start = crate::genome::G_OFFSET;
    let end = start + genome.n_genome as usize;
    f.write_all(&genome.g[start..end])?;
    Ok(())
}

/// Write the suffix-array PackedArray to `SA`.
pub fn write_sa(dir_out: &Path, sa: &PackedArray) -> Result<()> {
    let mut f = open_w(&dir_out.join("SA"))?;
    f.write_all(sa.as_bytes())?;
    Ok(())
}

/// Write the SAindex file: 1×`uint`(gSAindexNbases) + (gSAindexNbases+1)×`uint`
/// (`genomeSAindexStart`) + `SAi.charArray` bytes.
pub fn write_sai(
    dir_out: &Path,
    sa_index_nbases: u64,
    genome_sa_index_start: &[u64],
    sai: &PackedArray,
) -> Result<()> {
    let mut f = open_w(&dir_out.join("SAindex"))?;
    f.write_all(&sa_index_nbases.to_le_bytes())?;
    for v in genome_sa_index_start {
        f.write_all(&v.to_le_bytes())?;
    }
    f.write_all(sai.as_bytes())?;
    Ok(())
}

/// `genomeParametersWrite(fileName, P, errorOut, mapGen)`.
pub fn genome_parameters_write(
    file_name: &Path,
    p: &Parameters,
    genome: &Genome,
) -> Result<()> {
    let mut f = open_w(file_name)?;
    writeln!(f, "### {}", p.command_line_full)?;
    writeln!(f, "### GstrandBit {}", genome.g_strand_bit as i32)?;
    writeln!(f, "versionGenome\t{}", p.version_genome)?;
    writeln!(f, "genomeType\t{}", p.p_ge.g_type_string)?;

    write!(f, "genomeFastaFiles\t")?;
    for ff in &p.p_ge.g_fasta_files {
        write!(f, "{} ", ff)?;
    }
    writeln!(f)?;
    writeln!(f, "genomeSAindexNbases\t{}", p.p_ge.g_sa_index_nbases)?;
    writeln!(f, "genomeChrBinNbits\t{}", p.p_ge.g_chr_bin_nbits)?;
    writeln!(f, "genomeSAsparseD\t{}", p.p_ge.g_sa_sparse_d)?;

    writeln!(
        f,
        "genomeTransformType\t{}",
        p.p_ge.transform.ty_string
    )?;
    writeln!(f, "genomeTransformVCF\t{}", p.p_ge.transform.vcf_file)?;

    writeln!(f, "sjdbOverhang\t{}", genome.sjdb_overhang)?;
    write!(f, "sjdbFileChrStartEnd\t")?;
    for s in &p.p_ge.sjdb_file_chr_start_end {
        write!(f, "{} ", s)?;
    }
    writeln!(f)?;

    writeln!(f, "sjdbGTFfile\t{}", p.p_ge.sjdb_gtf_file)?;
    writeln!(f, "sjdbGTFchrPrefix\t{}", p.p_ge.sjdb_gtf_chr_prefix)?;
    writeln!(f, "sjdbGTFfeatureExon\t{}", p.p_ge.sjdb_gtf_feature_exon)?;
    writeln!(
        f,
        "sjdbGTFtagExonParentTranscript\t{}",
        p.p_ge.sjdb_gtf_tag_exon_parent_transcript
    )?;
    writeln!(
        f,
        "sjdbGTFtagExonParentGene\t{}",
        p.p_ge.sjdb_gtf_tag_exon_parent_gene
    )?;
    writeln!(f, "sjdbInsertSave\t{}", p.p_ge.sjdb_insert_save)?;

    // genomeFileSizes is space-separated after the key.
    write!(f, "genomeFileSizes\t{}", p.p_ge.g_file_sizes[0])?;
    for v in &p.p_ge.g_file_sizes[1..] {
        write!(f, " {}", v)?;
    }
    writeln!(f)?;

    Ok(())
}

fn open_w(path: &Path) -> Result<BufWriter<File>> {
    let f = File::create(path).with_context(|| format!("could not open {}", path.display()))?;
    Ok(BufWriter::new(f))
}
