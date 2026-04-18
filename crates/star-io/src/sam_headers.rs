//! 1:1 port of `samHeaders.cpp` for the SAM text output path.
//!
//! BAM header writing (via noodles BGZF) is deferred to M7. Transcriptome
//! header generation is deferred to M7 as well.

use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use star_core::STAR_VERSION;
use star_genome::Genome;
use star_params::parameters::Parameters;

/// Port of `samHeaders` (SAM-only subset). Returns the final `@HD` +
/// body header string (without the sort-order tag) and the
/// SO:coordinate variant for BAM-sorted output.
///
/// Writes the SAM header into `sam_out` if `P.outSAMbool`.
pub fn write_sam_headers(
    p: &Parameters,
    genome: &Genome,
    sam_out: &mut impl Write,
) -> anyhow::Result<(String, String)> {
    if p.out_sam_mode == "None" || p.out_sam_type.first().map(|s| s.as_str()) == Some("None") {
        return Ok((String::new(), String::new()));
    }

    let mut body = String::new();
    for ii in 0..genome.n_chr_real as usize {
        body.push_str(&format!(
            "@SQ\tSN:{}\tLN:{}\n",
            genome.chr_name[ii], genome.chr_length[ii]
        ));
    }

    let extra_path: PathBuf = format!("{}/extraReferences.txt", p.p_ge.g_dir).into();
    if let Ok(f) = fs::File::open(&extra_path) {
        let reader = BufReader::new(f);
        for line in reader.lines() {
            let line = line?;
            if !line.trim().is_empty() {
                body.push_str(&line);
                body.push('\n');
            }
        }
    }

    body.push_str(&format!(
        "@PG\tID:STAR\tPN:STAR\tVN:{}\tCL:{}\n",
        STAR_VERSION, p.command_line_full
    ));

    body.push_str(&format!("@CO\tuser command line: {}\n", p.command_line));

    let header_hd = "@HD\tVN:1.4".to_string();
    let sam_header = format!("{}\n{}", header_hd, body);
    let sam_header_sorted_coord = format!("{}\tSO:coordinate\n{}", header_hd, body);

    // Header is written when SAM is requested or when BAM output uses an
    // intermediate SAM file as scratch; the BAM path needs a valid SAM
    // header so that the SAM→BAM converter can extract @SQ entries.
    let needs_header = p.out_sam_type.first().map(|s| s.as_str()) == Some("SAM")
        || p.out_sam_type.iter().any(|s| s == "SAM")
        || p.out_bam_unsorted
        || p.out_bam_coord;
    if needs_header {
        sam_out.write_all(sam_header.as_bytes())?;
    }

    Ok((sam_header, sam_header_sorted_coord))
}
