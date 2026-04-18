//! Port of `source/Parameters.{h,cpp}` and `ParametersGenome.{h,cpp}`.
//!
//! Current scope: M2 — only fields required by `genomeGenerate` are wired to
//! the CLI parser. Other fields exist with default values for forward-compat
//! so higher layers can typecheck without NPEs.

use std::collections::HashSet;

use star_core::STAR_VERSION;

/// Port of `class ParametersGenome` (see `ParametersGenome.h`).
#[derive(Debug, Clone)]
pub struct ParametersGenome {
    /// `--genomeDir` / `gDir`
    pub g_dir: String,
    /// `--genomeLoad` / `gLoad`
    pub g_load: String,
    /// `gType` / `gTypeString`
    pub g_type: u32,
    pub g_type_string: String,
    /// `--genomeFastaFiles`
    pub g_fasta_files: Vec<String>,
    /// `--genomeChainFiles`
    pub g_chain_files: Vec<String>,
    /// Genome-transform sub-struct.
    pub transform: Transform,
    /// `--genomeSAindexNbases`
    pub g_sa_index_nbases: u64,
    /// `--genomeChrBinNbits`
    pub g_chr_bin_nbits: u64,
    /// `--genomeSAsparseD`
    pub g_sa_sparse_d: u64,
    /// `--genomeSuffixLengthMax`
    pub g_suffix_length_max: u64,
    /// `--genomeFileSizes` (populated on genomeGenerate).
    pub g_file_sizes: Vec<u64>,

    /// `--sjdbFileChrStartEnd`
    pub sjdb_file_chr_start_end: Vec<String>,
    /// `--sjdbGTFfile`
    pub sjdb_gtf_file: String,
    /// `--sjdbGTFchrPrefix`
    pub sjdb_gtf_chr_prefix: String,
    /// `--sjdbGTFfeatureExon`
    pub sjdb_gtf_feature_exon: String,
    /// `--sjdbGTFtagExonParentTranscript`
    pub sjdb_gtf_tag_exon_parent_transcript: String,
    /// `--sjdbGTFtagExonParentGene`
    pub sjdb_gtf_tag_exon_parent_gene: String,
    pub sjdb_gtf_tag_exon_parent_gene_name: Vec<String>,
    pub sjdb_gtf_tag_exon_parent_gene_type: Vec<String>,

    pub sjdb_insert_save: String,
    pub sjdb_overhang: u64,
    pub sjdb_overhang_par: i32,
    pub sjdb_score: i32,

    pub chr_set: ChrSet,
}

#[derive(Debug, Clone, Default)]
pub struct Transform {
    pub ty: i32,
    pub ty_string: String,
    pub vcf_file: String,
    pub output: Vec<String>,
    pub out_yes: bool,
    pub out_sam: bool,
    pub out_sj: bool,
    pub out_quant: bool,
}

#[derive(Debug, Clone, Default)]
pub struct ChrSet {
    pub mito_strings: Vec<String>,
    pub mito: HashSet<u64>,
}

/// Port of `Parameters::alignEndsType` nested struct
/// (Parameters.h:117-120).
///
/// `in` holds the textual option; `ext[iFrag][end]` flags whether that end
/// should be forced into end-to-end extension (as opposed to local clipping).
#[derive(Debug, Clone)]
pub struct AlignEndsType {
    pub in_: String,
    /// `ext[iFrag][end]` - end: 0 = 5', 1 = 3'.
    pub ext: [[bool; 2]; 2],
}

impl Default for AlignEndsType {
    fn default() -> Self {
        Self {
            in_: "Local".to_string(),
            ext: [[false; 2]; 2],
        }
    }
}

impl AlignEndsType {
    /// Port of the finalization block in `Parameters::inputParameters()`
    /// (Parameters.cpp:966-989).
    pub fn finalize(&mut self) -> anyhow::Result<()> {
        self.ext = [[false; 2]; 2];
        match self.in_.as_str() {
            "EndToEnd" => {
                self.ext = [[true; 2]; 2];
            }
            "Extend5pOfRead1" => self.ext[0][0] = true,
            "Extend5pOfReads12" => {
                self.ext[0][0] = true;
                self.ext[1][0] = true;
            }
            "Extend3pOfRead1" => self.ext[0][1] = true,
            "Local" => {}
            _ => anyhow::bail!(
                "EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --alignEndsType: {}",
                self.in_
            ),
        }
        Ok(())
    }
}

/// Port of `Parameters::alignEndsProtrude` (Parameters.h:122-126).
#[derive(Debug, Clone)]
pub struct AlignEndsProtrude {
    pub in_: Vec<String>,
    pub n_bases_max: i32,
    pub concordant_pair: bool,
}

impl Default for AlignEndsProtrude {
    fn default() -> Self {
        Self {
            in_: vec!["0".to_string(), "ConcordantPair".to_string()],
            n_bases_max: 0,
            concordant_pair: true,
        }
    }
}

impl AlignEndsProtrude {
    /// Port of Parameters.cpp:1084-1096.
    pub fn finalize(&mut self) -> anyhow::Result<()> {
        self.n_bases_max = self
            .in_
            .first()
            .map(|s| s.parse::<i32>())
            .transpose()?
            .unwrap_or(0);
        self.concordant_pair = false;
        if self.n_bases_max > 0 {
            let mode = self.in_.get(1).map(|s| s.as_str()).unwrap_or("");
            match mode {
                "ConcordantPair" => self.concordant_pair = true,
                "DiscordantPair" => self.concordant_pair = false,
                other => anyhow::bail!(
                    "EXITING because of fatal PARAMETERS error: unrecognized option in of --alignEndsProtrude={}",
                    other
                ),
            }
        }
        Ok(())
    }
}

/// Port of `Parameters::alignInsertionFlush` (Parameters.h:128-131).
#[derive(Debug, Clone)]
pub struct AlignInsertionFlush {
    pub in_: String,
    pub flush_right: bool,
}

impl Default for AlignInsertionFlush {
    fn default() -> Self {
        Self {
            in_: "None".to_string(),
            flush_right: false,
        }
    }
}

impl AlignInsertionFlush {
    /// Port of Parameters.cpp:1099-1108.
    pub fn finalize(&mut self) -> anyhow::Result<()> {
        match self.in_.as_str() {
            "None" => self.flush_right = false,
            "Right" => self.flush_right = true,
            other => anyhow::bail!(
                "EXITING because of fatal PARAMETERS error: unrecognized option in of --alignInsertionFlush={}",
                other
            ),
        }
        Ok(())
    }
}

impl Default for ParametersGenome {
    fn default() -> Self {
        Self {
            g_dir: "./GenomeDir/".to_string(),
            g_load: "NoSharedMemory".to_string(),
            g_type: 0,
            g_type_string: "Full".to_string(),
            g_fasta_files: Vec::new(),
            g_chain_files: Vec::new(),
            transform: Transform {
                ty_string: "None".to_string(),
                vcf_file: "-".to_string(),
                ..Default::default()
            },
            g_sa_index_nbases: 14,
            g_chr_bin_nbits: 18,
            g_sa_sparse_d: 1,
            g_suffix_length_max: u64::MAX,
            g_file_sizes: Vec::new(),
            sjdb_file_chr_start_end: vec!["-".to_string()],
            sjdb_gtf_file: "-".to_string(),
            sjdb_gtf_chr_prefix: "-".to_string(),
            sjdb_gtf_feature_exon: "exon".to_string(),
            sjdb_gtf_tag_exon_parent_transcript: "transcript_id".to_string(),
            sjdb_gtf_tag_exon_parent_gene: "gene_id".to_string(),
            sjdb_gtf_tag_exon_parent_gene_name: vec!["gene_name".to_string()],
            sjdb_gtf_tag_exon_parent_gene_type: vec![
                "gene_type".to_string(),
                "gene_biotype".to_string(),
            ],
            sjdb_insert_save: "Basic".to_string(),
            sjdb_overhang: 100,
            sjdb_overhang_par: 0,
            sjdb_score: 2,
            chr_set: ChrSet::default(),
        }
    }
}

/// Port of `class Parameters` (partial). Only fields needed so far; others
/// added as we grow the port. Field names match the C++ for grep-ability.
#[derive(Debug, Clone)]
pub struct Parameters {
    /// `commandLine` (joined args string).
    pub command_line: String,
    /// `commandLineFull` (written to genomeParameters.txt).
    pub command_line_full: String,
    /// `runMode`
    pub run_mode: String,
    pub run_thread_n: i32,
    pub run_dir_perm: u32,
    pub run_dir_perm_in: String,
    pub run_rng_seed: u64,

    /// `versionGenome`
    pub version_genome: String,

    /// `outFileNamePrefix`
    pub out_file_name_prefix: String,
    pub out_tmp_dir: String,
    pub out_log_file_name: String,

    pub p_ge: ParametersGenome,

    /// `limitGenomeGenerateRAM`
    pub limit_genome_generate_ram: u64,
    /// `limitSjdbInsertNsj`
    pub limit_sjdb_insert_nsj: u64,
    pub limit_io_buffer_size: u64,
    pub limit_out_sam_one_read_bytes: u64,
    pub limit_out_sj_collapsed: u64,
    pub limit_out_sj_one_read: u64,
    pub limit_ba_msort_ram: u64,
    pub limit_ba_msort_bytes: u64,
    pub limit_n_re_ad_ssj: u64,

    // ---- M3 alignment-time parameters (Parameters.{h,cpp}) ----
    /// `--alignIntronMax`
    pub align_intron_max: u64,
    pub align_intron_min: u64,
    pub align_mates_gap_max: u64,
    pub align_sj_overhang_min: u64,
    pub align_sj_stitch_mismatch_nmax: [i64; 4],
    pub align_sjdb_overhang_min: u64,
    pub align_sj_db_overhang_min: u64,
    pub align_splice_mismatch_nmax: u32,
    pub align_transcripts_per_read_nmax: u32,
    pub align_window_per_read_nmax: u32,
    pub align_transcripts_per_window_nmax: u32,
    pub align_ends_type: AlignEndsType,
    pub align_ends_protrude: AlignEndsProtrude,
    pub align_soft_clip_at_reference_ends: String,
    pub align_insertion_flush: AlignInsertionFlush,
    pub align_ends_protrude_nbases: u64,
    pub align_ends_protrude_conc_nonoverlap: bool,
    pub pe_overlap_nbases_min: i32,
    pub pe_overlap_mmp: f64,

    // Window / bin parameters (winBinNbits, winBinN, winFlankNbins, ...).
    pub win_bin_nbits: u64,
    pub win_bin_n: u64,
    pub win_bin_chr_nbits: u64,
    pub win_flank_nbins: u64,
    pub win_anchor_dist_nbins: u64,
    pub win_anchor_multimap_nmax: u64,
    pub win_read_coverage_relative_min: f64,
    pub win_read_coverage_basic_min: u64,

    // Output SAM/BAM related.
    pub out_sam_type: Vec<String>,
    /// Derived from `out_sam_type` by [`Parameters::finalize`]: true when
    /// `--outSAMtype SAM` is in effect. Mirrors `Parameters::outSAMbool`.
    pub out_sam_bool: bool,
    /// Derived: `--outSAMtype BAM Unsorted` in effect.
    pub out_bam_unsorted: bool,
    /// Derived: `--outSAMtype BAM SortedByCoordinate` in effect.
    pub out_bam_coord: bool,
    /// `--outBAMcompression` (1..9 or -1..=9 per STAR's conventions).
    pub out_bam_compression: i32,
    pub out_sam_mode: String,
    pub out_sam_strand_field: String,
    pub out_sam_attributes: Vec<String>,
    pub out_sam_unmapped: String,
    pub out_sam_order: String,
    pub out_sam_primary_flag: String,
    pub out_sam_read_id: String,
    pub out_sam_mapq_unique: u32,
    pub out_sam_flag_or: u32,
    pub out_sam_flag_and: u32,
    pub out_sam_heading: String,
    pub out_sam_multi_nmax: i64,
    pub out_sam_filter: Vec<String>,
    pub out_sam_tl_en: String,
    pub out_multimapper_order: String,
    pub out_multimapper_order_random: bool,

    // Read input parameters.
    pub read_files_in: Vec<String>,
    pub read_files_command: Vec<String>,
    pub read_files_type: Vec<String>,
    pub read_files_prefix: String,
    pub read_strand: String,
    pub read_map_number: i64,
    pub read_ends_length: i32,
    pub read_mates_length_same: bool,
    pub read_name_separator: Vec<String>,
    pub read_quality_score_base: u32,

    // Output filter / scoring.
    pub out_filter_type: String,
    pub out_filter_multimap_nmax: u32,
    pub out_filter_multimap_score_range: u32,
    pub out_filter_score_min: i32,
    pub out_filter_score_min_over_lread: f64,
    pub out_filter_match_nmin: i32,
    pub out_filter_match_nmin_over_lread: f64,
    pub out_filter_mismatch_nmax: i32,
    pub out_filter_mismatch_nover_lmax: f64,
    pub out_filter_mismatch_nover_rmax: f64,
    pub out_filter_intron_motifs: String,
    pub out_filter_intron_strands: String,

    // Score parameters.
    pub score_ins_open: i32,
    pub score_ins_base: i32,
    pub score_del_open: i32,
    pub score_del_base: i32,
    pub score_gap: i32,
    pub score_gap_gcag: i32,
    pub score_gap_atac: i32,
    pub score_gap_noncan: i32,
    pub score_genomic_length_log2_scale: f64,
    pub score_match: i32,
    pub score_mismatch: i32,
    pub score_stitch_sj_shift: i32,

    // Tracks whether user supplied sjdbOverhang on the CLI. Mirrors C++
    // `parArray.at(pGe.sjdbOverhang_par)->inputLevel > 0`.
    pub sjdb_overhang_user_set: bool,

    // ---- Seed parameters (Parameters.h:134-143) ----
    pub seed_multimap_nmax: u64,
    pub seed_search_lmax: u64,
    pub seed_per_read_nmax: u64,
    pub seed_per_window_nmax: u64,
    pub seed_none_loci_per_window: u64,
    pub seed_search_start_lmax: u64,
    pub seed_search_start_lmax_over_lread: f64,
    pub seed_split_min: u64,
    pub seed_map_min: u64,
    /// `maxNsplit` — hard cap on the number of good regions returned by
    /// `qualitySplit`. Mirrors Parameters.cpp:473 (value 10).
    pub max_n_split: u64,
    /// `readNends` — number of input streams (1 = SE, 2 = PE, 3 = PE + barcode).
    pub read_nends: u32,
    /// `readNmates` — number of mate sequences that contribute to alignment
    /// (excludes barcodes).
    pub read_nmates: u32,

    // ---- alignSplicedMateMapLmin / outFilterBySJoutStage / outSJ ----
    pub align_spliced_mate_map_lmin: u64,
    pub align_spliced_mate_map_lmin_over_lmate: f64,
    pub out_filter_by_sjout_stage: i32,

    // Novel SJs loaded from external filter — used by outFilterBySJoutStage=2 path.
    pub sj_novel_start: Vec<u64>,
    pub sj_novel_end: Vec<u64>,
    pub sj_novel_n: u64,

    // outSAMstrandField finalized value (type 0=none, 1=intronMotif).
    pub out_sam_strand_field_type: u32,

    // Chimeric minimum segment length (Parameters::pCh.segmentMin). 0 disables.
    // Kept as a top-level field for grep-ability; always mirrors `p_ch.segment_min`.
    pub p_ch_segment_min: u64,

    /// Port of `Parameters::pCh` (ParametersChimeric) — M6.
    pub p_ch: ParametersChimeric,

    /// Port of `Parameters::twoPass` (Parameters.h:258-266).
    pub two_pass: TwoPass,

    /// Port of `Parameters::sjdbInsert` (Parameters.h:269-274).
    pub sjdb_insert: SjdbInsert,

    /// Port of `Parameters::quant` (ParametersQuant.h + Parameters.cpp:891-947).
    pub quant: ParametersQuant,
}

/// Port of the anonymous `twoPass` struct inside `Parameters`.
///
/// Mirrors Parameters.cpp:267-269, 779-826. `mode` is the text value of
/// `--twopassMode` ("None" or "Basic"); `yes` is derived during
/// `finalize()`.
#[derive(Debug, Clone)]
pub struct TwoPass {
    pub yes: bool,
    pub pass2: bool,
    /// `--twopass1readsN`: int64 to accommodate the C++ default value of
    /// `-1` which signifies "map all reads in pass1".
    pub pass1_reads_n: i64,
    /// `true` if `--twopass1readsN` was explicitly set on the CLI. Used by
    /// `finalize` to diagnose the combination with `--twopassMode None`.
    pub pass1_reads_n_set: bool,
    /// Directory holding pass1 outputs; `<outFileNamePrefix>_STARpass1/`.
    pub dir: String,
    /// Path to the `SJ.out.tab` produced by pass1.
    pub pass1_sj_file: String,
    /// `"None"` (default) or `"Basic"`.
    pub mode: String,
}

impl Default for TwoPass {
    fn default() -> Self {
        Self {
            yes: false,
            pass2: false,
            pass1_reads_n: -1,
            pass1_reads_n_set: false,
            dir: String::new(),
            pass1_sj_file: String::new(),
            mode: "None".to_string(),
        }
    }
}

/// Port of the anonymous `sjdbInsert` struct inside `Parameters`.
#[derive(Debug, Clone, Default)]
pub struct SjdbInsert {
    pub yes: bool,
    pub pass1: bool,
    pub pass2: bool,
    /// `<outFileNamePrefix>_STARgenome/` — directory that holds the
    /// runtime genome produced when junctions are inserted on the fly.
    pub out_dir: String,
}

/// Port of `Parameters::quant.trSAM` (ParametersQuant.h lines 14-24 +
/// Parameters.cpp:891-935).
#[derive(Debug, Clone)]
pub struct QuantTrSAM {
    pub yes: bool,
    pub bam_yes: bool,
    pub bam_compression: i32,
    pub output: String,
    pub indel: bool,
    pub soft_clip: bool,
    pub single_end: bool,
}

impl Default for QuantTrSAM {
    fn default() -> Self {
        Self {
            yes: false,
            bam_yes: false,
            bam_compression: 1,
            output: "BanSingleEnd_BanIndels_ExtendSoftclip".to_string(),
            indel: false,
            soft_clip: false,
            single_end: false,
        }
    }
}

/// Port of `Parameters::quant.geCount` (GeneCounts sub-struct).
#[derive(Debug, Clone, Default)]
pub struct QuantGeCount {
    pub yes: bool,
    pub out_file: String,
}

/// Port of `Parameters::quant` (ParametersQuant.h + Parameters.cpp:891-947).
///
/// This port currently wires the two most commonly used sub-modes:
/// - `TranscriptomeSAM` (→ `Aligned.toTranscriptome.out.bam`)
/// - `GeneCounts`       (→ `ReadsPerGene.out.tab`)
///
/// The `geneFull*`, `gene`, and `sjdbGenome` sub-modes are not yet ported —
/// they will be added alongside their respective output paths.
#[derive(Debug, Clone, Default)]
pub struct ParametersQuant {
    /// `--quantMode` text (`-`, `GeneCounts`, `TranscriptomeSAM`, ...).
    pub mode: Vec<String>,
    /// Derived: any quant mode active (`mode[0] != "-"`).
    pub yes: bool,
    /// TranscriptomeSAM-related settings.
    pub tr_sam: QuantTrSAM,
    /// GeneCounts-related settings.
    pub ge_count: QuantGeCount,
}

/// Port of `ParametersChimeric` (ParametersChimeric.h + ParametersChimeric_initialize.cpp).
///
/// All `*in` text-valued CLI fields that map to bools (`out.type`, `filter.stringIn`)
/// are stored alongside the resolved bool; the resolved fields are derived
/// during [`ParametersChimeric::initialize`].
#[derive(Debug, Clone)]
pub struct ParametersChimeric {
    /// `--chimSegmentMin`. 0 disables chimeric detection entirely.
    pub segment_min: u64,
    /// `--chimJunctionOverhangMin`.
    pub junction_overhang_min: u64,
    /// `--chimSegmentReadGapMax`.
    pub segment_read_gap_max: u64,
    /// `--chimScoreMin`.
    pub score_min: i32,
    /// `--chimScoreDropMax`.
    pub score_drop_max: i32,
    /// `--chimScoreSeparation`.
    pub score_separation: i32,
    /// `--chimScoreJunctionNonGTAG`.
    pub score_junction_non_gtag: i32,
    /// `--chimMainSegmentMultNmax`.
    pub main_segment_mult_nmax: u64,

    /// `--chimMultimapScoreRange`.
    pub multimap_score_range: u64,
    /// `--chimMultimapNmax`. 0 selects the "Old" algorithm; >0 selects "Mult".
    pub multimap_nmax: u64,
    /// `--chimNonchimScoreDropMin`.
    pub nonchim_score_drop_min: u64,

    /// `--chimOutJunctionFormat`.
    pub out_junction_format: Vec<i32>,

    /// `--chimFilter` (text values: `None`, `banGenomicN`).
    pub filter_string_in: Vec<String>,
    /// Derived from `filter_string_in` — true if `banGenomicN` is present.
    pub filter_genomic_n: bool,

    /// `--chimOutType` (text values: `Junctions`, `SeparateSAMold`,
    /// `WithinBAM`, `HardClip`, `SoftClip`).
    pub out_type: Vec<String>,
    pub out_bam: bool,
    pub out_bam_hard_clip: bool,
    pub out_sam_old: bool,
    pub out_junctions: bool,
}

impl Default for ParametersChimeric {
    fn default() -> Self {
        // Mirrors parametersDefault; `initialize()` resolves derived bools.
        Self {
            segment_min: 0,
            junction_overhang_min: 20,
            segment_read_gap_max: 0,
            score_min: 0,
            score_drop_max: 20,
            score_separation: 10,
            score_junction_non_gtag: -1,
            main_segment_mult_nmax: 10,
            multimap_score_range: 1,
            multimap_nmax: 0,
            nonchim_score_drop_min: 20,
            out_junction_format: vec![0],
            filter_string_in: vec!["banGenomicN".to_string()],
            filter_genomic_n: false,
            out_type: vec!["Junctions".to_string()],
            out_bam: false,
            out_bam_hard_clip: true,
            out_sam_old: false,
            out_junctions: false,
        }
    }
}

impl ParametersChimeric {
    /// 1:1 port of `ParametersChimeric::initialize` (ParametersChimeric_initialize.cpp).
    ///
    /// Validates `--chimOutType` / `--chimFilter` and populates the derived
    /// `out_bam` / `out_sam_old` / `out_junctions` / `filter_genomic_n`
    /// fields. No files are opened here — output streams are created later
    /// by the caller if the relevant flag is set.
    pub fn initialize(&mut self, pe_overlap_nbases_min: u64) -> anyhow::Result<()> {
        self.out_bam = false;
        self.out_junctions = false;
        self.out_sam_old = false;
        self.out_bam_hard_clip = true;

        if self.segment_min == 0 {
            return Ok(());
        }

        for t in &self.out_type {
            match t.as_str() {
                "WithinBAM" => self.out_bam = true,
                "SeparateSAMold" => self.out_sam_old = true,
                "Junctions" => self.out_junctions = true,
                "HardClip" => self.out_bam_hard_clip = true,
                "SoftClip" => self.out_bam_hard_clip = false,
                other => anyhow::bail!(
                    "EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --chimOutType: {}\n\
                     SOLUTION: re-run STAR with --chimOutType Junctions , SeparateSAMold , WithinBAM , HardClip",
                    other
                ),
            }
        }

        if self.multimap_nmax > 0 && self.out_sam_old {
            anyhow::bail!(
                "EXITING because of fatal PARAMETERS error: --chimMultimapNmax > 0 (new chimeric detection) presently only works with --chimOutType Junctions/WithinBAM\n\
                 SOLUTION: re-run with --chimOutType Junctions/WithinBAM"
            );
        }

        if pe_overlap_nbases_min > 0
            && self.multimap_nmax == 0
            && (self.out_junctions || self.out_sam_old)
        {
            anyhow::bail!(
                "EXITING because of fatal PARAMETERS error: --chimMultimapNmax 0 (default old chimeric detection) and --peOverlapNbasesMin > 0 (merging overlapping mates) presently only works with --chimOutType WithinBAM\n\
                 SOLUTION: re-run with --chimOutType WithinBAM"
            );
        }

        self.filter_genomic_n = false;
        for f in &self.filter_string_in {
            match f.as_str() {
                "banGenomicN" => self.filter_genomic_n = true,
                "None" => {}
                other => anyhow::bail!(
                    "EXITING because of fatal PARAMETERS error: unrecognized value of --chimFilter={}\n\
                     SOLUTION: use allowed values: banGenomicN || None",
                    other
                ),
            }
        }
        Ok(())
    }
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            command_line: String::new(),
            command_line_full: String::new(),
            run_mode: "alignReads".to_string(),
            run_thread_n: 1,
            run_dir_perm: 0o700,
            run_dir_perm_in: "User_RWX".to_string(),
            run_rng_seed: 777,
            version_genome: "2.7.4a".to_string(),
            out_file_name_prefix: "./".to_string(),
            out_tmp_dir: "-".to_string(),
            out_log_file_name: "./Log.out".to_string(),
            p_ge: ParametersGenome::default(),
            limit_genome_generate_ram: 31_000_000_000,
            limit_sjdb_insert_nsj: 1_000_000,
            limit_io_buffer_size: 150_000_000,
            limit_out_sam_one_read_bytes: 100_000,
            limit_out_sj_collapsed: 1_000_000,
            limit_out_sj_one_read: 1_000,
            limit_ba_msort_ram: 0,
            limit_ba_msort_bytes: 0,
            limit_n_re_ad_ssj: 10_000_000,

            // Alignment / window defaults mirror parametersDefault.
            align_intron_max: 0,
            align_intron_min: 21,
            align_mates_gap_max: 0,
            align_sj_overhang_min: 5,
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            align_sjdb_overhang_min: 3,
            align_sj_db_overhang_min: 3,
            align_splice_mismatch_nmax: 0,
            align_transcripts_per_read_nmax: 10_000,
            align_window_per_read_nmax: 10_000,
            align_transcripts_per_window_nmax: 100,
            align_ends_type: AlignEndsType::default(),
            align_ends_protrude: AlignEndsProtrude::default(),
            align_soft_clip_at_reference_ends: "Yes".to_string(),
            align_insertion_flush: AlignInsertionFlush::default(),
            align_ends_protrude_nbases: 0,
            align_ends_protrude_conc_nonoverlap: true,
            pe_overlap_nbases_min: 0,
            pe_overlap_mmp: 0.01,

            win_bin_nbits: 16,
            win_bin_n: 0,
            win_bin_chr_nbits: 0,
            win_flank_nbins: 4,
            win_anchor_dist_nbins: 9,
            win_anchor_multimap_nmax: 50,
            win_read_coverage_relative_min: 0.5,
            win_read_coverage_basic_min: 0,

            out_sam_type: vec!["SAM".to_string()],
            out_sam_bool: true,
            out_bam_unsorted: false,
            out_bam_coord: false,
            out_bam_compression: 1,
            out_sam_mode: "Full".to_string(),
            out_sam_strand_field: "None".to_string(),
            out_sam_attributes: vec!["Standard".to_string()],
            out_sam_unmapped: "None".to_string(),
            out_sam_order: "Paired".to_string(),
            out_sam_primary_flag: "OneBestScore".to_string(),
            out_sam_read_id: "Standard".to_string(),
            out_sam_mapq_unique: 255,
            out_sam_flag_or: 0,
            out_sam_flag_and: 65535,
            out_sam_heading: "-".to_string(),
            out_sam_multi_nmax: -1,
            out_sam_filter: vec!["-".to_string()],
            out_sam_tl_en: "-".to_string(),
            out_multimapper_order: "Old_2.4".to_string(),
            out_multimapper_order_random: false,

            read_files_in: vec!["Read1".to_string(), "Read2".to_string()],
            read_files_command: vec!["-".to_string()],
            read_files_type: vec!["Fastx".to_string()],
            read_files_prefix: "-".to_string(),
            read_strand: "Unstranded".to_string(),
            read_map_number: -1,
            read_ends_length: -1,
            read_mates_length_same: false,
            read_name_separator: vec!["/".to_string()],
            read_quality_score_base: 33,

            out_filter_type: "Normal".to_string(),
            out_filter_multimap_nmax: 10,
            out_filter_multimap_score_range: 1,
            out_filter_score_min: 0,
            out_filter_score_min_over_lread: 0.66,
            out_filter_match_nmin: 0,
            out_filter_match_nmin_over_lread: 0.66,
            out_filter_mismatch_nmax: 10,
            out_filter_mismatch_nover_lmax: 0.3,
            out_filter_mismatch_nover_rmax: 1.0,
            out_filter_intron_motifs: "None".to_string(),
            out_filter_intron_strands: "RemoveInconsistentStrands".to_string(),

            score_ins_open: -2,
            score_ins_base: -2,
            score_del_open: -2,
            score_del_base: -2,
            score_gap: 0,
            score_gap_gcag: -4,
            score_gap_atac: -8,
            score_gap_noncan: -8,
            score_genomic_length_log2_scale: -0.25,
            score_match: 1,
            score_mismatch: -1,
            // C++ parametersDefault / Parameters.cpp: the factory default is 1,
            // not 200. With 200 the move-left phase of stitchAlignToTranscript's
            // motif scan explores ~60 extra positions and finds alignments that
            // C++ STAR never considers (observed on Ath SRR6281254 read
            // SRR6281254.sra.12786635 where Rust produced a spurious 5-exon
            // 47kb-intron candidate beating the correct 3-exon alignment).
            score_stitch_sj_shift: 1,

            sjdb_overhang_user_set: false,

            seed_multimap_nmax: 10_000,
            seed_search_lmax: 0,
            seed_per_read_nmax: 1_000,
            seed_per_window_nmax: 50,
            seed_none_loci_per_window: 10,
            seed_search_start_lmax: 50,
            seed_search_start_lmax_over_lread: 1.0,
            seed_split_min: 12,
            seed_map_min: 5,
            max_n_split: 10,
            read_nends: 1,
            read_nmates: 1,

            align_spliced_mate_map_lmin: 0,
            align_spliced_mate_map_lmin_over_lmate: 0.66,
            out_filter_by_sjout_stage: 0,

            sj_novel_start: Vec::new(),
            sj_novel_end: Vec::new(),
            sj_novel_n: 0,

            out_sam_strand_field_type: 0,
            p_ch_segment_min: 0,
            p_ch: ParametersChimeric::default(),
            two_pass: TwoPass::default(),
            sjdb_insert: SjdbInsert::default(),
            quant: ParametersQuant {
                mode: vec!["-".to_string()],
                yes: false,
                tr_sam: QuantTrSAM::default(),
                ge_count: QuantGeCount::default(),
            },
        }
    }
}

impl Parameters {
    pub fn new() -> Self {
        Self::default()
    }

    /// Alias matching the C++ builder pattern; right now identical to [`new`].
    pub fn new_with_defaults() -> Self {
        Self::default()
    }

    /// Simplified command-line parser for milestones M1..M2. Accepts the
    /// multi-valued `--<key> v1 v2 ...` style used by STAR. Unknown keys cause
    /// an error per the reference behavior.
    pub fn parse_cli(&mut self, args: &[String]) -> anyhow::Result<()> {
        self.command_line = args.join(" ");
        self.command_line_full = self.command_line.clone();
        let mut i = 1; // skip argv[0]
        while i < args.len() {
            let key = &args[i];
            if !key.starts_with("--") {
                anyhow::bail!(
                    "EXITING because of fatal PARAMETERS error: unrecognized parameter '{key}'"
                );
            }
            let key = &key[2..];
            let mut values = Vec::new();
            i += 1;
            while i < args.len() && !args[i].starts_with("--") {
                values.push(args[i].clone());
                i += 1;
            }
            self.set_param(key, &values)?;
        }
        self.finalize()?;
        Ok(())
    }

    /// Port of finalization logic from `Parameters::inputParameters()`:
    /// resolves textual option strings into typed fields.
    pub fn finalize(&mut self) -> anyhow::Result<()> {
        self.align_ends_type.finalize()?;
        self.align_ends_protrude.finalize()?;
        self.align_insertion_flush.finalize()?;

        // Port of Parameters_readFilesInit.cpp:152-166 — for default FASTX
        // input (`readFilesTypeN==1`), the number of mates equals the number
        // of input files (1 = SE, 2 = PE). `readNmates` is the same (Solo
        // may later reduce it by 1, but we don't support Solo in M4).
        if self.run_mode == "alignReads" && !self.read_files_in.is_empty() {
            let n = self.read_files_in.len() as u32;
            if (1..=2).contains(&n) {
                self.read_nends = n;
                self.read_nmates = n;
            }
        }

        self.finalize_two_pass()?;
        self.finalize_sjdb_insert()?;
        self.finalize_out_sam_type()?;

        // Keep the legacy mirror in sync with the structured field.
        self.p_ch_segment_min = self.p_ch.segment_min;
        // ParametersChimeric::initialize (M6): validates --chimOutType /
        // --chimFilter and populates derived bools. Must be called AFTER
        // peOverlap is parsed so we can pass `pe_overlap_nbases_min`.
        self.p_ch
            .initialize(self.pe_overlap_nbases_min.max(0) as u64)?;

        self.finalize_quant()?;
        Ok(())
    }

    /// Port of Parameters.cpp:608-683 — `--outSAMtype` dispatch. Sets
    /// `out_sam_bool`, `out_bam_unsorted`, `out_bam_coord`.
    fn finalize_out_sam_type(&mut self) -> anyhow::Result<()> {
        self.out_sam_bool = false;
        self.out_bam_unsorted = false;
        self.out_bam_coord = false;
        if self.run_mode != "alignReads" || self.out_sam_mode == "None" {
            return Ok(());
        }
        if self.out_sam_type.is_empty() {
            return Ok(());
        }
        match self.out_sam_type[0].as_str() {
            "BAM" => {
                if self.out_sam_type.len() < 2 {
                    anyhow::bail!(
                        "EXITING because of fatal PARAMETER error: missing BAM option\n\
                         SOLUTION: re-run STAR with one of the allowed values of --outSAMtype BAM Unsorted OR SortedByCoordinate OR both"
                    );
                }
                for kind in &self.out_sam_type[1..] {
                    match kind.as_str() {
                        "Unsorted" => self.out_bam_unsorted = true,
                        "SortedByCoordinate" => self.out_bam_coord = true,
                        other => anyhow::bail!(
                            "EXITING because of fatal input ERROR: unknown value for --outSAMtype: {other}"
                        ),
                    }
                }
            }
            "SAM" => {
                if self.out_sam_type.len() > 1 {
                    anyhow::bail!(
                        "EXITING because of fatal PARAMETER error: --outSAMtype SAM cannot be combined with {}",
                        self.out_sam_type[1]
                    );
                }
                self.out_sam_bool = true;
            }
            "None" => {}
            other => anyhow::bail!(
                "EXITING because of fatal input ERROR: unknown value for the first word of --outSAMtype: {other}"
            ),
        }
        Ok(())
    }

    /// Port of Parameters.cpp:891-947 — `quant` finalisation.
    fn finalize_quant(&mut self) -> anyhow::Result<()> {
        self.quant.yes = false;
        self.quant.ge_count.yes = false;
        self.quant.tr_sam.yes = false;
        self.quant.tr_sam.bam_yes = false;
        self.quant.tr_sam.indel = false;
        self.quant.tr_sam.soft_clip = false;
        self.quant.tr_sam.single_end = false;

        if self.quant.mode.is_empty() || self.quant.mode[0] == "-" {
            return Ok(());
        }
        self.quant.yes = true;
        for m in self.quant.mode.clone() {
            match m.as_str() {
                "-" => {}
                "TranscriptomeSAM" => {
                    self.quant.tr_sam.yes = true;
                    if self.quant.tr_sam.bam_compression > -2 {
                        self.quant.tr_sam.bam_yes = true;
                    }
                    match self.quant.tr_sam.output.as_str() {
                        "BanSingleEnd_BanIndels_ExtendSoftclip" => {
                            self.quant.tr_sam.indel = false;
                            self.quant.tr_sam.soft_clip = false;
                            self.quant.tr_sam.single_end = false;
                        }
                        "BanSingleEnd" => {
                            self.quant.tr_sam.indel = true;
                            self.quant.tr_sam.soft_clip = true;
                            self.quant.tr_sam.single_end = false;
                        }
                        "BanSingleEnd_ExtendSoftclip" => {
                            self.quant.tr_sam.indel = true;
                            self.quant.tr_sam.soft_clip = false;
                            self.quant.tr_sam.single_end = false;
                        }
                        other => anyhow::bail!(
                            "EXITING because of fatal INPUT error: unknown value of --quantTranscriptomeSAMoutput={other}"
                        ),
                    }
                }
                "GeneCounts" => {
                    self.quant.ge_count.yes = true;
                    self.quant.ge_count.out_file =
                        format!("{}ReadsPerGene.out.tab", self.out_file_name_prefix);
                }
                other => anyhow::bail!(
                    "EXITING because of fatal INPUT error: unrecognized option in --quantMode={other}\n\
                     SOLUTION: use one of the allowed values of --quantMode : TranscriptomeSAM or GeneCounts or - ."
                ),
            }
        }
        Ok(())
    }

    /// Port of Parameters.cpp:779-826 — `twoPass` finalisation.
    fn finalize_two_pass(&mut self) -> anyhow::Result<()> {
        if self.two_pass.pass1_reads_n_set && self.two_pass.mode == "None" {
            anyhow::bail!(
                "EXITING because of fatal PARAMETERS error: --twopass1readsN is defined, but --twoPassMode is not defined\n\
                 SOLUTION: to activate the 2-pass mode, use --twopassMode Basic"
            );
        }

        self.two_pass.yes = false;
        self.two_pass.pass2 = false;

        if self.two_pass.mode != "None" {
            if self.run_mode != "alignReads" {
                anyhow::bail!(
                    "EXITING because of fatal PARAMETERS error: 2-pass mapping is not compatible with runMode={}\n\
                     SOLUTION: remove --twopassMode option",
                    self.run_mode
                );
            }
            if self.two_pass.mode != "Basic" {
                anyhow::bail!(
                    "EXITING because of fatal PARAMETERS error: unrecognized value of --twopassMode={}\n\
                     SOLUTION: for the 2-pass mode, use allowed values --twopassMode: Basic",
                    self.two_pass.mode
                );
            }
            if self.two_pass.pass1_reads_n == 0 {
                anyhow::bail!(
                    "EXITING because of fatal PARAMETERS error: --twopass1readsN = 0 in the 2-pass mode\n\
                     SOLUTION: for the 2-pass mode, specify --twopass1readsN > 0. Use a very large number or -1 to map all reads in the 1st pass."
                );
            }
            if self.p_ge.g_load != "NoSharedMemory" {
                anyhow::bail!(
                    "EXITING because of fatal PARAMETERS error: 2-pass mode is not compatible with --genomeLoad={}\n\
                     SOLUTION: re-run STAR with --genomeLoad NoSharedMemory ; this is the only option compatible with --twopassMode Basic .",
                    self.p_ge.g_load
                );
            }
            self.two_pass.yes = true;
            self.two_pass.dir = format!("{}_STARpass1/", self.out_file_name_prefix);
            std::fs::create_dir_all(&self.two_pass.dir).map_err(|e| {
                anyhow::anyhow!(
                    "EXITING because of fatal ERROR: could not make pass1 directory: {}: {e}",
                    self.two_pass.dir
                )
            })?;
        }
        Ok(())
    }

    /// Port of Parameters.cpp:999-1035 — `sjdbInsert` finalisation.
    fn finalize_sjdb_insert(&mut self) -> anyhow::Result<()> {
        self.sjdb_insert.pass1 = false;
        self.sjdb_insert.pass2 = false;
        self.sjdb_insert.yes = false;

        let has_sjdb_file = !self.p_ge.sjdb_file_chr_start_end.is_empty()
            && self.p_ge.sjdb_file_chr_start_end[0] != "-";
        let has_sjdb_gtf = self.p_ge.sjdb_gtf_file != "-";
        if has_sjdb_file || has_sjdb_gtf {
            self.sjdb_insert.pass1 = true;
            self.sjdb_insert.yes = true;
        }
        if self.two_pass.yes {
            self.sjdb_insert.pass2 = true;
            self.sjdb_insert.yes = true;
        }

        if self.p_ge.g_load != "NoSharedMemory" && self.sjdb_insert.yes {
            anyhow::bail!(
                "EXITING because of fatal PARAMETERS error: on the fly junction insertion and 2-pass mapping cannot be used with shared memory genome\n\
                 SOLUTION: run STAR with --genomeLoad NoSharedMemory to avoid using shared memory"
            );
        }

        if self.run_mode == "alignReads" && self.sjdb_insert.yes {
            if self.p_ge.sjdb_overhang == 0 {
                anyhow::bail!(
                    "EXITING because of fatal PARAMETERS error: pGe.sjdbOverhang <=0 while junctions are inserted on the fly with --sjdbFileChrStartEnd or/and --sjdbGTFfile\n\
                     SOLUTION: specify pGe.sjdbOverhang>0, ideally readmateLength-1"
                );
            }
            self.sjdb_insert.out_dir = format!("{}_STARgenome/", self.out_file_name_prefix);
            std::fs::create_dir_all(&self.sjdb_insert.out_dir).map_err(|e| {
                anyhow::anyhow!(
                    "EXITING because of fatal ERROR: could not make run-time genome directory: {}: {e}",
                    self.sjdb_insert.out_dir
                )
            })?;
        }
        Ok(())
    }

    fn set_param(&mut self, key: &str, values: &[String]) -> anyhow::Result<()> {
        let single = || -> anyhow::Result<&str> {
            values
                .first()
                .map(|s| s.as_str())
                .ok_or_else(|| anyhow::anyhow!("--{key}: missing value"))
        };
        match key {
            "runMode" => self.run_mode = single()?.to_string(),
            "runThreadN" => self.run_thread_n = single()?.parse()?,
            "runRNGseed" => self.run_rng_seed = single()?.parse()?,
            "runDirPerm" => self.run_dir_perm_in = single()?.to_string(),
            "outFileNamePrefix" => self.out_file_name_prefix = single()?.to_string(),
            "outTmpDir" => self.out_tmp_dir = single()?.to_string(),
            "genomeDir" => self.p_ge.g_dir = single()?.to_string(),
            "genomeLoad" => self.p_ge.g_load = single()?.to_string(),
            "genomeFastaFiles" => {
                self.p_ge.g_fasta_files = values.to_vec();
            }
            "genomeChainFiles" => self.p_ge.g_chain_files = values.to_vec(),
            "genomeSAindexNbases" => self.p_ge.g_sa_index_nbases = single()?.parse()?,
            "genomeChrBinNbits" => self.p_ge.g_chr_bin_nbits = single()?.parse()?,
            "genomeSAsparseD" => self.p_ge.g_sa_sparse_d = single()?.parse()?,
            "genomeSuffixLengthMax" => self.p_ge.g_suffix_length_max = single()?.parse()?,
            "sjdbFileChrStartEnd" => {
                self.p_ge.sjdb_file_chr_start_end = values.to_vec();
            }
            "sjdbGTFfile" => self.p_ge.sjdb_gtf_file = single()?.to_string(),
            "sjdbGTFchrPrefix" => self.p_ge.sjdb_gtf_chr_prefix = single()?.to_string(),
            "sjdbGTFfeatureExon" => self.p_ge.sjdb_gtf_feature_exon = single()?.to_string(),
            "sjdbGTFtagExonParentTranscript" => {
                self.p_ge.sjdb_gtf_tag_exon_parent_transcript = single()?.to_string()
            }
            "sjdbGTFtagExonParentGene" => {
                self.p_ge.sjdb_gtf_tag_exon_parent_gene = single()?.to_string()
            }
            "sjdbOverhang" => {
                self.p_ge.sjdb_overhang = single()?.parse()?;
                self.p_ge.sjdb_overhang_par = 1;
                self.sjdb_overhang_user_set = true;
            }
            "sjdbScore" => self.p_ge.sjdb_score = single()?.parse()?,
            "sjdbInsertSave" => self.p_ge.sjdb_insert_save = single()?.to_string(),
            "limitGenomeGenerateRAM" => self.limit_genome_generate_ram = single()?.parse()?,
            "limitSjdbInsertNsj" => self.limit_sjdb_insert_nsj = single()?.parse()?,
            "versionGenome" => self.version_genome = single()?.to_string(),
            "twopassMode" => self.two_pass.mode = single()?.to_string(),
            "twopass1readsN" => {
                self.two_pass.pass1_reads_n = single()?.parse()?;
                self.two_pass.pass1_reads_n_set = true;
            }

            // --- M6 chimeric detection (ParametersChimeric) ---
            "chimOutType" => self.p_ch.out_type = values.to_vec(),
            "chimSegmentMin" => {
                self.p_ch.segment_min = single()?.parse()?;
                self.p_ch_segment_min = self.p_ch.segment_min;
            }
            "chimScoreMin" => self.p_ch.score_min = single()?.parse()?,
            "chimScoreDropMax" => self.p_ch.score_drop_max = single()?.parse()?,
            "chimScoreSeparation" => self.p_ch.score_separation = single()?.parse()?,
            "chimScoreJunctionNonGTAG" => self.p_ch.score_junction_non_gtag = single()?.parse()?,
            "chimJunctionOverhangMin" => self.p_ch.junction_overhang_min = single()?.parse()?,
            "chimSegmentReadGapMax" => self.p_ch.segment_read_gap_max = single()?.parse()?,
            "chimFilter" => self.p_ch.filter_string_in = values.to_vec(),
            "chimMainSegmentMultNmax" => self.p_ch.main_segment_mult_nmax = single()?.parse()?,
            "chimMultimapNmax" => self.p_ch.multimap_nmax = single()?.parse()?,
            "chimMultimapScoreRange" => self.p_ch.multimap_score_range = single()?.parse()?,
            "chimNonchimScoreDropMin" => self.p_ch.nonchim_score_drop_min = single()?.parse()?,
            "chimOutJunctionFormat" => {
                self.p_ch.out_junction_format =
                    values.iter().filter_map(|s| s.parse().ok()).collect();
            }

            // --- M7 quantification (ParametersQuant) ---
            "quantMode" => self.quant.mode = values.to_vec(),
            "quantTranscriptomeBAMcompression" => {
                self.quant.tr_sam.bam_compression = single()?.parse()?
            }
            "quantTranscriptomeSAMoutput" => self.quant.tr_sam.output = single()?.to_string(),
            // Accept-and-ignore for M2 (these don't affect genome generation).
            "genomeTransformType" => self.p_ge.transform.ty_string = single()?.to_string(),
            "genomeTransformVCF" => self.p_ge.transform.vcf_file = single()?.to_string(),
            "parametersFiles" | "sysShell" => {
                // deferred to later milestones
            }
            "readFilesIn" => self.read_files_in = values.to_vec(),
            "readFilesCommand" => self.read_files_command = values.to_vec(),
            "readFilesType" => self.read_files_type = values.to_vec(),
            "readFilesPrefix" => self.read_files_prefix = single()?.to_string(),
            "readMapNumber" => self.read_map_number = single()?.parse()?,

            // --- M3 alignment params ---
            "alignIntronMin" => self.align_intron_min = single()?.parse()?,
            "alignIntronMax" => self.align_intron_max = single()?.parse()?,
            "alignMatesGapMax" => self.align_mates_gap_max = single()?.parse()?,
            "alignSJoverhangMin" => self.align_sj_overhang_min = single()?.parse()?,
            "alignSJDBoverhangMin" => {
                self.align_sjdb_overhang_min = single()?.parse()?;
                self.align_sj_db_overhang_min = self.align_sjdb_overhang_min;
            }
            "alignSplicedMateMapLmin" => self.align_spliced_mate_map_lmin = single()?.parse()?,
            "alignSplicedMateMapLminOverLmate" => {
                self.align_spliced_mate_map_lmin_over_lmate = single()?.parse()?
            }
            "seedMultimapNmax" => self.seed_multimap_nmax = single()?.parse()?,
            "seedSearchLmax" => self.seed_search_lmax = single()?.parse()?,
            "seedPerReadNmax" => self.seed_per_read_nmax = single()?.parse()?,
            "seedPerWindowNmax" => self.seed_per_window_nmax = single()?.parse()?,
            "seedNoneLociPerWindow" => self.seed_none_loci_per_window = single()?.parse()?,
            "seedSearchStartLmax" => self.seed_search_start_lmax = single()?.parse()?,
            "seedSearchStartLmaxOverLread" => {
                self.seed_search_start_lmax_over_lread = single()?.parse()?
            }
            "seedSplitMin" => self.seed_split_min = single()?.parse()?,
            "seedMapMin" => self.seed_map_min = single()?.parse()?,
            "outFilterBySJoutStage" => self.out_filter_by_sjout_stage = single()?.parse()?,
            "alignWindowsPerReadNmax" => self.align_window_per_read_nmax = single()?.parse()?,
            "alignTranscriptsPerWindowNmax" => {
                self.align_transcripts_per_window_nmax = single()?.parse()?
            }
            "alignTranscriptsPerReadNmax" => {
                self.align_transcripts_per_read_nmax = single()?.parse()?
            }
            "alignEndsType" => self.align_ends_type.in_ = single()?.to_string(),
            "alignEndsProtrude" => self.align_ends_protrude.in_ = values.to_vec(),
            "alignInsertionFlush" => self.align_insertion_flush.in_ = single()?.to_string(),
            "alignSoftClipAtReferenceEnds" => {
                self.align_soft_clip_at_reference_ends = single()?.to_string()
            }
            "winBinNbits" => self.win_bin_nbits = single()?.parse()?,
            "winAnchorDistNbins" => self.win_anchor_dist_nbins = single()?.parse()?,
            "winFlankNbins" => self.win_flank_nbins = single()?.parse()?,
            "winAnchorMultimapNmax" => self.win_anchor_multimap_nmax = single()?.parse()?,
            "winReadCoverageRelativeMin" => {
                self.win_read_coverage_relative_min = single()?.parse()?
            }
            "winReadCoverageBasicMin" => self.win_read_coverage_basic_min = single()?.parse()?,

            "outSAMtype" => self.out_sam_type = values.to_vec(),
            "outBAMcompression" => self.out_bam_compression = single()?.parse()?,
            "outSAMmode" => self.out_sam_mode = single()?.to_string(),
            "outSAMstrandField" => self.out_sam_strand_field = single()?.to_string(),
            "outSAMattributes" => self.out_sam_attributes = values.to_vec(),
            "outSAMunmapped" => self.out_sam_unmapped = values.to_vec().join(" "),
            "outSAMorder" => self.out_sam_order = single()?.to_string(),
            "outSAMprimaryFlag" => self.out_sam_primary_flag = single()?.to_string(),
            "outSAMreadID" => self.out_sam_read_id = single()?.to_string(),
            "outSAMmapqUnique" => self.out_sam_mapq_unique = single()?.parse()?,
            "outSAMmultNmax" => self.out_sam_multi_nmax = single()?.parse()?,
            "outSAMheaderHD" | "outSAMheaderPG" | "outSAMheaderCommentFile" => {}
            "outMultimapperOrder" => self.out_multimapper_order = single()?.to_string(),

            "outFilterType" => self.out_filter_type = single()?.to_string(),
            "outFilterMultimapNmax" => self.out_filter_multimap_nmax = single()?.parse()?,
            "outFilterMultimapScoreRange" => {
                self.out_filter_multimap_score_range = single()?.parse()?
            }
            "outFilterScoreMin" => self.out_filter_score_min = single()?.parse()?,
            "outFilterScoreMinOverLread" => {
                self.out_filter_score_min_over_lread = single()?.parse()?
            }
            "outFilterMatchNmin" => self.out_filter_match_nmin = single()?.parse()?,
            "outFilterMatchNminOverLread" => {
                self.out_filter_match_nmin_over_lread = single()?.parse()?
            }
            "outFilterMismatchNmax" => self.out_filter_mismatch_nmax = single()?.parse()?,
            "outFilterMismatchNoverLmax" => {
                self.out_filter_mismatch_nover_lmax = single()?.parse()?
            }
            "outFilterMismatchNoverReadLmax" => {
                self.out_filter_mismatch_nover_rmax = single()?.parse()?
            }
            "outFilterIntronMotifs" => self.out_filter_intron_motifs = single()?.to_string(),
            "outFilterIntronStrands" => self.out_filter_intron_strands = single()?.to_string(),

            "scoreInsOpen" => self.score_ins_open = single()?.parse()?,
            "scoreInsBase" => self.score_ins_base = single()?.parse()?,
            "scoreDelOpen" => self.score_del_open = single()?.parse()?,
            "scoreDelBase" => self.score_del_base = single()?.parse()?,
            "scoreGap" => self.score_gap = single()?.parse()?,
            "scoreGapGCAG" => self.score_gap_gcag = single()?.parse()?,
            "scoreGapATAC" => self.score_gap_atac = single()?.parse()?,
            "scoreGapNoncan" => self.score_gap_noncan = single()?.parse()?,
            "scoreGenomicLengthLog2scale" => {
                self.score_genomic_length_log2_scale = single()?.parse()?
            }
            "scoreStitchSJshift" => self.score_stitch_sj_shift = single()?.parse()?,
            _ => {
                // For M1/M2 we accept unknown keys silently so the CLI layer
                // doesn't fail on align-path parameters we haven't wired yet.
                // TODO: restore `exitWithError` behavior once full Parameters
                // port lands (M3).
            }
        }
        Ok(())
    }
}

/// Build-time compilation info; equivalent to `COMPILATION_TIME_PLACE`.
pub fn compilation_time_place() -> &'static str {
    env!("STAR_RS_COMPILATION_TIME_PLACE")
}

pub fn star_version() -> &'static str {
    STAR_VERSION
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_parses_genome_generate() {
        let args: Vec<String> = [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            "/tmp/gd",
            "--genomeFastaFiles",
            "/tmp/a.fa",
            "/tmp/b.fa",
            "--genomeSAindexNbases",
            "11",
            "--runThreadN",
            "4",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let mut p = Parameters::new();
        p.parse_cli(&args).unwrap();
        assert_eq!(p.run_mode, "genomeGenerate");
        assert_eq!(p.p_ge.g_dir, "/tmp/gd");
        assert_eq!(p.p_ge.g_fasta_files, vec!["/tmp/a.fa", "/tmp/b.fa"]);
        assert_eq!(p.p_ge.g_sa_index_nbases, 11);
        assert_eq!(p.run_thread_n, 4);
    }
}
