//! 1:1 port of `source/STAR.cpp` (top-level `main()`).
//!
//! During M1, only the CLI dispatch, `usage()`, and the very first phase of
//! parameter parsing are ported. Subsequent milestones flesh out the branches.

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use star_params::parameters::{compilation_time_place, star_version};
use star_params::PARAMETERS_DEFAULT;

/// `void usage(int usageType)` from `STAR.cpp`.
fn usage(usage_type: i32) -> ! {
    let stdout = std::io::stdout();
    let mut w = stdout.lock();
    let _ = writeln!(
        w,
        "Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq"
    );
    let _ = writeln!(
        w,
        "Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2022\n"
    );
    let _ = writeln!(w, "STAR version={}", star_version());
    let _ = writeln!(
        w,
        "STAR compilation time,server,dir={}",
        compilation_time_place()
    );
    let _ = writeln!(w, "For more details see:");
    let _ = writeln!(w, "<https://github.com/alexdobin/STAR>");
    let _ = writeln!(
        w,
        "<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>"
    );
    if usage_type == 0 {
        let _ = writeln!(w, "\nTo list all parameters, run STAR --help");
    } else if usage_type == 1 {
        let _ = w.write_all(PARAMETERS_DEFAULT);
    }
    drop(w);
    std::process::exit(0);
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let argc = args.len();

    if argc == 1 {
        usage(0);
    }
    if argc == 2 && (args[1] == "-h" || args[1] == "--help") {
        usage(1);
    }

    let mut p = star_params::parameters::Parameters::new_with_defaults();
    if let Err(e) = p.parse_cli(&args) {
        eprintln!("EXITING: fatal error parsing parameters: {}", e);
        return ExitCode::from(1);
    }
    p.command_line_full = args[..].join(" ");

    match p.run_mode.as_str() {
        "genomeGenerate" => match run_genome_generate(&mut p) {
            Ok(_) => ExitCode::from(0),
            Err(e) => {
                eprintln!("EXITING because of FATAL ERROR in genomeGenerate: {e:#}");
                ExitCode::from(102)
            }
        },
        "alignReads" => match run_align_reads(&mut p) {
            Ok(_) => ExitCode::from(0),
            Err(e) => {
                eprintln!("EXITING because of FATAL ERROR in alignReads: {e:#}");
                ExitCode::from(103)
            }
        },
        other => {
            eprintln!("EXITING: unsupported --runMode {other}");
            ExitCode::from(101)
        }
    }
}

fn run_genome_generate(p: &mut star_params::parameters::Parameters) -> anyhow::Result<()> {
    let mut genome = star_genome::genome_generate(p)?;
    if p.sjdb_insert.yes {
        p.sjdb_insert.out_dir = p.p_ge.g_dir.clone();
        if genome.chr_bin.is_empty() {
            genome.chr_bin_fill();
        }

        let mut sjdb_loci = star_sjdb::SjdbLoci::new();
        let mut gtf = star_sjdb::Gtf::new(&mut genome, p, &p.p_ge.g_dir, &mut sjdb_loci)?;
        gtf.transcript_gene_sj(&genome, p, &p.p_ge.g_dir, &mut sjdb_loci)?;
        sjdb_loci.load_from_files(&p.p_ge.sjdb_file_chr_start_end)?;

        let genome1 = genome.clone();
        star_sjdb::insert_junctions::sjdb_insert_junctions(
            p,
            &mut genome,
            &genome1,
            &mut sjdb_loci,
        )?;

        p.p_ge.g_file_sizes = vec![genome.n_genome, genome.n_sa_byte];
        let dir = Path::new(&p.p_ge.g_dir);
        star_genome::io::write_genome_sequence(dir, &genome)?;
        star_genome::io::write_sa(dir, &genome.sa)?;
        star_genome::io::write_sai(
            dir,
            p.p_ge.g_sa_index_nbases,
            &genome.genome_sa_index_start,
            &genome.sai,
        )?;
        star_genome::io::genome_parameters_write(&dir.join("genomeParameters.txt"), p, &genome)?;
    }
    Ok(())
}

/// Top-level for `--runMode alignReads`. Implements the STAR.cpp
/// 2-pass driver: optional pass-1 (with restricted parameters) → sjdb
/// insertion → main mapping pass.
fn run_align_reads(p: &mut star_params::parameters::Parameters) -> anyhow::Result<()> {
    use std::time::{SystemTime, UNIX_EPOCH};

    let now = || {
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs() as i64)
            .unwrap_or(0)
    };
    let t_start = now();

    let mut genome = star_genome::Genome::new();
    star_genome::genome_load(p, &mut genome)?;

    // ----- 2-pass: pass 1 -----
    if p.two_pass.yes {
        run_pass1(p, &genome)?;
        // Wire up for pass 2: the main alignReads run becomes pass 2.
        p.two_pass.pass2 = true;
        p.two_pass.pass1_sj_file = format!("{}SJ.out.tab", p.two_pass.dir);
        p.sjdb_insert.pass2 = true;
        p.sjdb_insert.yes = true;
    }

    // ----- sjdb insertion (annotation + pass1 SJs) -----
    if p.sjdb_insert.yes {
        let genome1 = genome.clone();
        let mut sjdb_loci = star_sjdb::SjdbLoci::new();
        star_sjdb::insert_junctions::sjdb_insert_junctions(
            p,
            &mut genome,
            &genome1,
            &mut sjdb_loci,
        )?;
    }

    // ----- main mapping pass -----
    let t_start_map = now();
    let n_threads = p.run_thread_n.max(1) as usize;
    let prefix_owned = p.out_file_name_prefix.clone();
    let (n_reads, mut stats, chunk_out_sj) = run_mapping_pass(p, &genome, &prefix_owned, true)?;
    let t_finish = now();

    stats.time_start = t_start;
    stats.time_start_map = t_start_map;
    stats.time_finish_map = t_finish;
    stats.time_finish = t_finish;

    let prefix = p.out_file_name_prefix.clone();
    let sj_filter = star_sjdb::output_sj::SjFilter::default();
    let sj_path = format!("{prefix}SJ.out.tab");
    let mut final_sj = chunk_out_sj;
    star_sjdb::output_sj::output_sj(
        std::slice::from_mut(&mut final_sj),
        &sj_filter,
        &genome.chr_name,
        &genome.chr_start,
        &sj_path,
    )?;

    let log_final_path = format!("{prefix}Log.final.out");
    let log_file = std::fs::File::create(&log_final_path)?;
    let mut log_out = std::io::BufWriter::new(log_file);
    stats.report_final(&mut log_out)?;
    log_out.flush()?;

    eprintln!(
        "star-rs alignReads: processed {} reads -> {}Aligned.out.sam ({} SJs, {} threads)",
        n_reads,
        prefix,
        final_sj.n(),
        n_threads,
    );
    Ok(())
}

/// Port of `twoPassRunPass1.cpp`. Runs the mapping pipeline with a
/// restricted `Parameters` into `p.two_pass.dir`, producing a pass-1
/// `SJ.out.tab` to be consumed by the main pass.
fn run_pass1(
    p: &mut star_params::parameters::Parameters,
    genome: &star_genome::Genome,
) -> anyhow::Result<()> {
    // Clone P and disable output that is not needed for pass 1.
    let mut p1 = p.clone();
    p1.out_sam_type[0] = "None".to_string();
    p1.p_ch_segment_min = 0;
    // TODO: p1.quant.* = false once quant is wired in M7.
    p1.out_filter_by_sjout_stage = 0;
    p1.out_file_name_prefix = p.two_pass.dir.clone();
    if p.two_pass.pass1_reads_n > 0 {
        p1.read_map_number = std::cmp::min(p.two_pass.pass1_reads_n, p.read_map_number);
    }

    // Run mapping pass with SAM output suppressed (we don't need the SAM
    // for pass 1; only SJ.out.tab is consumed by the 2nd pass).
    let prefix1 = p1.out_file_name_prefix.clone();
    let (_n_reads, _stats, chunk_out_sj) = run_mapping_pass(&mut p1, genome, &prefix1, false)?;

    // Write pass1 SJ.out.tab to p.two_pass.dir/SJ.out.tab.
    let sj_filter = star_sjdb::output_sj::SjFilter::default();
    let sj_path = format!("{}SJ.out.tab", p.two_pass.dir);
    let mut sj = chunk_out_sj;
    star_sjdb::output_sj::output_sj(
        std::slice::from_mut(&mut sj),
        &sj_filter,
        &genome.chr_name,
        &genome.chr_start,
        &sj_path,
    )?;
    Ok(())
}

/// Runs one alignment pass (`ReadAlignChunk::process_chunks` +
/// SAM/BAM output), returning `(n_reads, stats, out_sj)` from the chunk.
///
/// When `emit_sam = false`, the SAM file is still opened (to keep
/// downstream code simple) but redirected to `/dev/null`. The
/// `Chimeric.out.junction` file (if `--chimOutType Junctions` is
/// enabled) is opened only when `emit_chim` is true.
fn run_mapping_pass(
    p: &mut star_params::parameters::Parameters,
    genome: &star_genome::Genome,
    out_prefix: &str,
    emit_sam: bool,
) -> anyhow::Result<(u64, star_stats::stats::Stats, star_sjdb::out_sj::OutSJ)> {
    use std::fs::File;
    use std::io::{BufReader, BufWriter};

    let sam_attrs = star_params::sam_attributes::SamAttributes::resolve(p)?;
    let sam_path = format!("{out_prefix}Aligned.out.sam");
    let bam_path = format!("{out_prefix}Aligned.out.bam");
    // Decide whether to keep the SAM on disk (`--outSAMtype SAM` or nothing
    // set — default), or to transiently produce SAM and then
    // convert/sort to BAM (`--outSAMtype BAM Unsorted` or
    // `--outSAMtype BAM SortedByCoordinate`). Phase-2a always goes via
    // the SAM intermediate.
    let want_sam = emit_sam && (p.out_sam_bool || (!p.out_bam_unsorted && !p.out_bam_coord));
    let want_bam_unsorted = emit_sam && p.out_bam_unsorted;
    let want_bam_coord = emit_sam && p.out_bam_coord;
    let want_bam = want_bam_unsorted || want_bam_coord;
    // When only BAM is requested, we still write SAM to a scratch path
    // (`<prefix>Aligned.out.sam.tmp`) as the intermediate, then convert
    // to BAM at the end and delete the intermediate.
    let sam_scratch_path = format!("{out_prefix}Aligned.out.sam.tmp");
    let sam_file: Box<dyn Write + Send> = if emit_sam {
        std::fs::create_dir_all(
            std::path::Path::new(out_prefix)
                .parent()
                .unwrap_or(std::path::Path::new(".")),
        )
        .ok();
        if want_sam {
            Box::new(BufWriter::new(File::create(&sam_path)?))
        } else if want_bam {
            Box::new(BufWriter::new(File::create(&sam_scratch_path)?))
        } else {
            Box::new(std::io::sink())
        }
    } else {
        Box::new(std::io::sink())
    };
    let mut sam_out = sam_file;

    if emit_sam {
        star_io::sam_headers::write_sam_headers(p, genome, &mut sam_out)?;
    }

    let mut readers: Vec<BufReader<File>> = Vec::new();
    for path in &p.read_files_in {
        let f = File::open(path)?;
        readers.push(BufReader::new(f));
    }
    if readers.is_empty() {
        anyhow::bail!("EXITING: --readFilesIn is empty");
    }

    let chim_enabled = emit_sam && p.p_ch.segment_min > 0 && p.p_ch.out_junctions;

    // Per-read chim hook: run chimericDetection + chimericDetectionOldOutput.
    let chim_hook = |pp: &star_params::parameters::Parameters,
                     gg: &star_genome::Genome,
                     ra: &mut star_align::read_align::ReadAlign,
                     buf: &mut Vec<u8>|
     -> anyhow::Result<()> {
        if pp.p_ch.segment_min == 0 {
            return Ok(());
        }
        star_chimeric::detection::chimeric_detection(ra, pp, gg)?;
        if ra.chim_record && pp.p_ch.out_junctions {
            star_chimeric::detection::chimeric_detection_old_output(ra, pp, gg, buf)?;
        }
        Ok(())
    };

    // M7.2 Transcriptome load: only if --quantMode GeneCounts (or
    // TranscriptomeSAM in M7.4+).
    let quant_enabled = emit_sam && p.quant.ge_count.yes;
    let transcriptome = if emit_sam && p.quant.yes {
        star_quant::transcriptome::Transcriptome::load(p)?
    } else {
        None
    };
    let n_ge = transcriptome.as_ref().map(|t| t.n_ge).unwrap_or(0);

    // Quant hook (runs per-read after chim hook, before output_alignments).
    // Mirrors `ReadAlign::alignedAnnotation` → `geneCountsAddAlign` in
    // `STAR/source/ReadAlign_outputAlignments.cpp:298-325`.
    let ex_g_opt = transcriptome.as_ref().map(|t| &t.ex_g);
    let quant_init = || star_quant::quant::Quantifications::new(n_ge);
    let quant_hook = |pp: &star_params::parameters::Parameters,
                      _gg: &star_genome::Genome,
                      ra: &mut star_align::read_align::ReadAlign,
                      q: &mut star_quant::quant::Quantifications|
     -> anyhow::Result<()> {
        if !pp.quant.ge_count.yes {
            return Ok(());
        }
        let Some(ex_g) = ex_g_opt else { return Ok(()) };
        // Equivalent of `if (unmapType<0)` gate in outputAlignments.
        if ra.unmap_type >= 0 {
            return Ok(());
        }
        // Build &[&Transcript] from tr_mult_array (the set of accepted alignments).
        let trs: Vec<&star_align::Transcript> = ra.tr_mult_array.iter().collect();
        let mut gene1 = [-1i32; 3];
        star_quant::quant::gene_counts_add_align(q, ex_g, &trs, &mut gene1);
        Ok(())
    };

    let n_threads = p.run_thread_n.max(1) as usize;
    let (n_reads, stats, chunk_out_sj, chim_bufs, quant_chunks) = if n_threads > 1 {
        let all_reads = star_align::threads::load_all_reads(p, &mut readers[..])?;
        star_align::threads::run_align_multithread(
            p,
            genome,
            all_reads,
            &sam_attrs.order,
            &mut sam_out,
            p.run_rng_seed as u32,
            n_threads,
            chim_hook,
            quant_init,
            quant_hook,
        )?
    } else {
        let mut chunk = star_align::chunk::ReadAlignChunk::new(
            0,
            p.run_rng_seed as u32,
            sam_attrs.order.clone(),
        );
        let mut q = quant_init();
        let n = chunk.process_chunks(
            p,
            genome,
            &mut readers[..],
            &mut sam_out,
            |ra, buf| chim_hook(p, genome, ra, buf),
            |ra| quant_hook(p, genome, ra, &mut q),
        )?;
        (
            n,
            chunk.stats_ra,
            chunk.chunk_out_sj,
            vec![chunk.chunk_out_chim_junction],
            vec![q],
        )
    };
    sam_out.flush()?;
    drop(sam_out);

    // M7.4 (minimum viable): if user requested BAM, convert the SAM
    // intermediate into Aligned.out.bam via noodles. Byte-exactness with
    // C++ STAR's custom alignBAM encoder is not guaranteed yet.
    if want_bam_unsorted {
        let compression_level = Some(p.out_bam_compression.clamp(0, 9) as u32);
        star_bam::sam_to_bam::convert_sam_file_to_bam(
            std::path::Path::new(&sam_scratch_path),
            std::path::Path::new(&bam_path),
            compression_level,
        )?;
        let _ = std::fs::remove_file(&sam_scratch_path);
    } else if want_bam_coord {
        let compression_level = Some(p.out_bam_compression.clamp(0, 9) as u32);
        let sorted_bam_path = format!("{out_prefix}Aligned.sortedByCoord.out.bam");
        let tmp_root = if p.out_bam_sort_tmp_dir.is_empty() {
            format!("{out_prefix}_STARtmp_BAMsort")
        } else {
            p.out_bam_sort_tmp_dir.clone()
        };
        std::fs::create_dir_all(&tmp_root)?;

        // Build a per-pass BamOutput, feed the SAM intermediate
        // through it, run the sort, clean up.
        let n_bins = p.out_bam_coord_nbins.max(2);
        let mut bam_out = star_bam::bam_output::BamOutput::new(
            0,
            std::path::Path::new(&tmp_root),
            n_bins,
            p.chunk_out_bam_size_bytes,
        )?;
        star_bam::bam_sort_feed::feed_sam_into_bam_output(
            std::path::Path::new(&sam_scratch_path),
            &mut bam_out,
        )?;
        bam_out.coord_flush()?;

        // Extract the SAM header bytes to hand to the final writer.
        let sam_header_bytes = std::fs::read(&sam_scratch_path)?;
        let mut hdr_end = 0usize;
        for line in sam_header_bytes.split_inclusive(|&b| b == b'\n') {
            if !line.starts_with(b"@") { break; }
            hdr_end += line.len();
        }
        let sam_header = &sam_header_bytes[..hdr_end];

        star_bam::bam_sort_coord::sort_bam_by_coordinate(
            std::path::Path::new(&tmp_root),
            /* n_threads */ 1,
            n_bins,
            &[bam_out.bin_total_n().to_vec()],
            &[bam_out.bin_total_bytes().to_vec()],
            sam_header,
            std::path::Path::new(&sorted_bam_path),
            compression_level,
            p.limit_bam_sort_ram,
        )?;
        let _ = std::fs::remove_file(&sam_scratch_path);
        let _ = std::fs::remove_dir_all(&tmp_root);
    }

    if chim_enabled {
        let chim_path = format!("{out_prefix}Chimeric.out.junction");
        let chim_file = File::create(&chim_path)?;
        let mut chim_out = BufWriter::new(chim_file);
        // Column header is emitted only when `multimapNmax>0`; the legacy
        // Old path omits it (matches ParametersChimeric_initialize.cpp:48-71).
        for buf in &chim_bufs {
            chim_out.write_all(buf)?;
        }
        chim_out.flush()?;
    }

    // M7.3: merge per-chunk Quantifications and write ReadsPerGene.out.tab.
    if quant_enabled {
        if let Some(tr) = transcriptome.as_ref() {
            let mut merged = star_quant::quant::Quantifications::new(tr.n_ge);
            for q in &quant_chunks {
                merged.add_quants(q);
            }
            let path = p.quant.ge_count.out_file.clone();
            write_reads_per_gene(&path, tr, &merged, &stats)?;
        }
    }

    Ok((n_reads, stats, chunk_out_sj))
}

/// Port of `Transcriptome::quantsOutput` (Transcriptome.cpp:156-190).
/// Writes `ReadsPerGene.out.tab` with 4 header rows (N_unmapped,
/// N_multimapping, N_noFeature, N_ambiguous) followed by per-gene rows.
fn write_reads_per_gene(
    path: &str,
    tr: &star_quant::transcriptome::Transcriptome,
    q: &star_quant::quant::Quantifications,
    stats: &star_stats::stats::Stats,
) -> anyhow::Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    // N_unmapped = unmappedMismatch + unmappedShort + unmappedOther + unmappedMulti
    let n_unmapped = stats.unmapped_mismatch
        + stats.unmapped_short
        + stats.unmapped_other
        + stats.unmapped_multi;
    let n_type = q.gene_counts.n_type;

    write!(w, "N_unmapped")?;
    for _ in 0..n_type {
        write!(w, "\t{}", n_unmapped)?;
    }
    writeln!(w)?;

    write!(w, "N_multimapping")?;
    for _ in 0..n_type {
        write!(w, "\t{}", q.gene_counts.c_multi)?;
    }
    writeln!(w)?;

    write!(w, "N_noFeature")?;
    for it in 0..n_type {
        write!(w, "\t{}", q.gene_counts.c_none[it])?;
    }
    writeln!(w)?;

    write!(w, "N_ambiguous")?;
    for it in 0..n_type {
        write!(w, "\t{}", q.gene_counts.c_ambig[it])?;
    }
    writeln!(w)?;

    for ig in 0..tr.n_ge as usize {
        write!(w, "{}", tr.ge_id[ig])?;
        for it in 0..n_type {
            write!(w, "\t{}", q.gene_counts.g_count[it][ig])?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}
