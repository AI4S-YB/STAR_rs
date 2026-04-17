//! 1:1 port of `genomeScanFastaFiles.cpp`.
//!
//! Two-pass scan: first with `flag_run=false` to compute chromosome metadata
//! and the total padded genome length `N`; second with `flag_run=true` to
//! actually copy the numeric codes into the genome buffer.

use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{anyhow, Context, Result};

use star_core::seq::convert_nucleotides_to_numbers_remove_controls;

use crate::genome::Genome;

/// `uint genomeScanFastaFiles(Parameters& P, char* G, bool flagRun, Genome&)`.
///
/// When `flag_run == false`, computes chromosome starts/lengths/names and
/// returns the padded total length `N`. `g` may be empty in this mode.
/// When `flag_run == true`, fills the numeric genome buffer `g` starting at
/// offset 0 and returns the final padded `N`.
pub fn genome_scan_fasta_files(
    fasta_files: &[String],
    genome_chr_bin_nbases: u64,
    g: &mut [u8],
    flag_run: bool,
    gen: &mut Genome,
    mut log_main: impl FnMut(&str),
) -> Result<u64> {
    let mut n: u64 = 0;
    // Port: "previous chr records exist" path — re-run after an insert.
    if !flag_run && !gen.chr_length.is_empty() {
        gen.chr_start.pop();
        n = *gen.chr_start.last().unwrap_or(&0) + *gen.chr_length.last().unwrap_or(&0);
        gen.chr_length.pop();
    }

    for (ii, path) in fasta_files.iter().enumerate() {
        let _ = ii;
        let f = File::open(path)
            .with_context(|| format!("EXITING because of INPUT ERROR: could not open genomeFastaFile: {path}"))?;
        let mut rdr = BufReader::new(f);

        // Peek first byte — must be '>'.
        let buf = rdr.fill_buf()?;
        if buf.is_empty() {
            return Err(anyhow!("could not read from genomeFastaFile: {path}"));
        }
        let cc = buf[0];
        if cc != b'>' {
            return Err(anyhow!(
                "EXITING because of INPUT ERROR: the file format of the genomeFastaFile: {} is not fasta: the first character is '{}' ({}), not '>'",
                path, cc as char, cc
            ));
        }

        let mut line = String::new();
        loop {
            line.clear();
            let n_read = rdr.read_line(&mut line)?;
            if n_read == 0 {
                break;
            }
            // Strip trailing newline(s) to match getline semantics.
            while matches!(line.as_bytes().last(), Some(b'\n') | Some(b'\r')) {
                line.pop();
            }
            let bytes = line.as_bytes();
            if bytes.first().copied() == Some(b'>') {
                if !flag_run {
                    // Port: lineInStream.ignore(1,' '); lineInStream >> chrName1;
                    // i.e. drop '>' then take first whitespace-delimited token.
                    let rest = &bytes[1..];
                    let tok_end = rest
                        .iter()
                        .position(|b| b.is_ascii_whitespace())
                        .unwrap_or(rest.len());
                    let chr_name = String::from_utf8_lossy(&rest[..tok_end]).into_owned();
                    gen.chr_name.push(chr_name);
                }

                if !flag_run && !gen.chr_start.is_empty() {
                    gen.chr_length.push(n - gen.chr_start.last().copied().unwrap_or(0));
                }

                if n > 0 {
                    n = ((n + 1) / genome_chr_bin_nbases + 1) * genome_chr_bin_nbases;
                }

                if !flag_run {
                    gen.chr_start.push(n);
                    let msg = format!(
                        "{} : chr # {}  \"{}\" chrStart: {}\n",
                        path,
                        gen.chr_start.len() - 1,
                        gen.chr_name.last().unwrap(),
                        n
                    );
                    log_main(&msg);
                }
            } else {
                if flag_run {
                    // C++: N += convertNucleotidesToNumbersRemoveControls(line, G+N, len)
                    let start = n as usize;
                    let end = start + bytes.len();
                    let slice = &mut g[start..end];
                    let written = convert_nucleotides_to_numbers_remove_controls(bytes, slice);
                    n += written;
                } else {
                    // Count non-control chars (>=32).
                    for &b in bytes {
                        if (b as i32) >= 32 {
                            n += 1;
                        }
                    }
                }
            }
        }
    }

    if !flag_run {
        let last_start = gen.chr_start.last().copied().unwrap_or(0);
        gen.chr_length.push(n - last_start);
    }

    n = ((n + 1) / genome_chr_bin_nbases + 1) * genome_chr_bin_nbases;

    if !flag_run {
        gen.n_chr_real = gen.chr_start.len() as u64;
        gen.chr_start.push(n); // size-of-genome sentinel
        log_main("Chromosome sequence lengths: \n");
        for ii in 0..gen.n_chr_real as usize {
            gen.chr_name_index
                .insert(gen.chr_name[ii].clone(), ii as u64);
            log_main(&format!("{}\t{}\n", gen.chr_name[ii], gen.chr_length[ii]));
        }
    }

    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    fn tempfile(dir: &tempdir::PathHelper, name: &str, contents: &str) -> PathBuf {
        let p = dir.0.join(name);
        let mut f = File::create(&p).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        p
    }

    mod tempdir {
        use std::path::PathBuf;
        pub struct PathHelper(pub PathBuf);
        impl PathHelper {
            pub fn new() -> Self {
                let p = std::env::temp_dir().join(format!("star-rs-{}", std::process::id()));
                let _ = std::fs::create_dir_all(&p);
                Self(p)
            }
        }
        impl Drop for PathHelper {
            fn drop(&mut self) {
                let _ = std::fs::remove_dir_all(&self.0);
            }
        }
    }

    #[test]
    fn scan_two_chromosomes() {
        let dir = tempdir::PathHelper::new();
        let fa = tempfile(
            &dir,
            "small.fa",
            ">chr1 description\nACGTACGT\nNNNN\n>chr2\nAAAA\n",
        );
        let chr_bin_nbases = 1u64 << 6; // 64
        let mut g = Genome::default();
        // Scan pass 1: sizes only.
        let n = genome_scan_fasta_files(
            &[fa.display().to_string()],
            chr_bin_nbases,
            &mut [],
            false,
            &mut g,
            |_| {},
        )
        .unwrap();
        assert_eq!(g.n_chr_real, 2);
        assert_eq!(g.chr_name, vec!["chr1".to_string(), "chr2".to_string()]);
        assert_eq!(g.chr_length[0], 12); // ACGTACGTNNNN
        assert_eq!(g.chr_length[1], 4);

        // Scan pass 2: fill sequence.
        let mut buf = vec![0u8; n as usize + 16];
        let mut g2 = Genome::default();
        let _ = genome_scan_fasta_files(
            &[fa.display().to_string()],
            chr_bin_nbases,
            &mut buf,
            false,
            &mut g2,
            |_| {},
        )
        .unwrap();
        let n2 = genome_scan_fasta_files(
            &[fa.display().to_string()],
            chr_bin_nbases,
            &mut buf,
            true,
            &mut g2,
            |_| {},
        )
        .unwrap();
        assert_eq!(n2, n);
    }
}
