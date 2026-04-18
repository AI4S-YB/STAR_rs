//! 1:1 port of `sjdbInsertJunctions.cpp`.
//!
//! Orchestrates: load sjdb → prepare → buildIndex. Reused for:
//! - runtime junction insertion via `--sjdbFileChrStartEnd` / `--sjdbGTFfile`
//! - `--twopassMode Basic` pass-2 setup (loads pass-1 SJ.out.tab).

use std::fs::File;
use std::io::BufReader;

use star_genome::genome::Genome;
use star_params::parameters::Parameters;

use crate::build_index::sjdb_build_index;
use crate::prepare::sjdb_prepare;
use crate::sjdb_class::SjdbLoci;

/// Port of `sjdbInsertJunctions` (sjdbInsertJunctions.cpp:11-102).
///
/// - `map_gen`: genome to be extended in-place.
/// - `map_gen1`: pre-insertion genome snapshot (must be cloned by caller
///   before `map_gen` was modified). Used to re-index old junctions.
/// - `sjdb_loci`: accumulated junction loci (GTF/file + pass1 SJ).
///
/// Side effects:
/// - `map_gen.sjdb_*` and SA / SAi are rewritten.
/// - Writes `sjdbInfo.txt` + `sjdbList.out.tab` into `p.sjdb_insert.out_dir`.
/// - Updates `p.win_bin_n` per Parameters.cpp:101.
pub fn sjdb_insert_junctions(
    p: &mut Parameters,
    map_gen: &mut Genome,
    map_gen1: &Genome,
    sjdb_loci: &mut SjdbLoci,
) -> anyhow::Result<()> {
    // 1. If genome already had saved sjdb and sjdbLoci is empty, load from
    //    `sjdbList.out.tab` in the genome dir.
    if map_gen.sjdb_n > 0 && sjdb_loci.is_empty() {
        let path = format!("{}/sjdbList.out.tab", p.p_ge.g_dir);
        let f = File::open(&path).map_err(|e| {
            anyhow::anyhow!(
                "EXITING because of fatal INPUT error: could not open {path}: {e}\n\
                 SOLUTION: re-generate the genome in pGe.gDir={}",
                p.p_ge.g_dir
            )
        })?;
        sjdb_loci.load_from_stream(BufReader::new(f))?;
        sjdb_loci.priority.resize(sjdb_loci.chr.len(), 30);
    }

    // 2. At pass-2, load the novel junctions produced by pass-1.
    if p.two_pass.pass2 {
        let path = p.two_pass.pass1_sj_file.clone();
        sjdb_loci.load_from_path_with_priority(&path, 0).map_err(|e| {
            anyhow::anyhow!(
                "FATAL INPUT error, could not open input file with junctions from the 1st pass={path}: {e}"
            )
        })?;
    } else if p.run_mode != "genomeGenerate" {
        sjdb_loci.load_from_files(&p.p_ge.sjdb_file_chr_start_end)?;
        // GTF loading is deferred (M5.9 covers external files; GTF is M7 quant).
    }

    // 3. Prepare the merged sjdb and build new junction sequences.
    let n_genome_real = map_gen.chr_start[map_gen.n_chr_real as usize];
    let mut gsj = sjdb_prepare(sjdb_loci, p, n_genome_real, &p.sjdb_insert.out_dir, map_gen)?;

    if map_gen.sjdb_n > p.limit_sjdb_insert_nsj {
        anyhow::bail!(
            "Fatal LIMIT error: the number of junctions to be inserted on the fly ={} is larger than the limitSjdbInsertNsj={}\n\
             SOLUTION: re-run with at least --limitSjdbInsertNsj {}",
            map_gen.sjdb_n, p.limit_sjdb_insert_nsj, map_gen.sjdb_n
        );
    }

    // 4. Build the new SA / SAi.
    sjdb_build_index(p, &mut gsj, map_gen, map_gen1)?;

    // 5. Update win_bin_n to reflect the new n_genome.
    let win_bin_nbits = p.win_bin_nbits;
    p.win_bin_n = map_gen.n_genome / (1u64 << win_bin_nbits) + 1;

    Ok(())
}
