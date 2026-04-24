//! 1:1 port of `GTF.cpp` + `GTF_transcriptGeneSJ.cpp`.
//!
//! The `Gtf` loader parses exon records from `--sjdbGTFfile`, collects
//! transcript / gene metadata, writes the transcriptome sidecar tables used by
//! `GeneCounts`, and appends collapsed GTF-derived junctions into `SjdbLoci`.

use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use star_core::compression::open_maybe_compressed;
use star_genome::genome::Genome;
use star_params::parameters::Parameters;

use crate::sjdb_class::SjdbLoci;

pub const EX_T: usize = 0;
pub const EX_S: usize = 1;
pub const EX_E: usize = 2;
pub const EX_G: usize = 3;
pub const EX_L: usize = 4;

const EXGE_EX_START: usize = 0;
const EXGE_EX_END: usize = 1;
const EXGE_EX_STRAND: usize = 2;
const EXGE_GE_ID: usize = 3;
const EXGE_TR_ID: usize = 4;
const EXGE_LEN: usize = 5;

const EXTR_TR_START: usize = 0;
const EXTR_TR_END: usize = 1;
const EXTR_TR_ID: usize = 2;
const EXTR_EX_START: usize = 3;
const EXTR_EX_END: usize = 4;
const EXTR_GE_ID: usize = 5;
const EXTR_LEN: usize = 6;

#[derive(Debug, Default, Clone)]
pub struct Gtf {
    pub gtf_yes: bool,
    pub exon_loci: Vec<[u64; EX_L]>,
    pub transcript_strand: Vec<u8>,
    pub transcript_id: Vec<String>,
    pub gene_id: Vec<String>,
    pub gene_attr: Vec<[String; 2]>,
}

impl Gtf {
    pub fn new(
        genome: &mut Genome,
        p: &Parameters,
        dir_out: impl AsRef<Path>,
        _sjdb_loci: &mut SjdbLoci,
    ) -> Result<Self> {
        if genome.sjdb_overhang == 0 || p.p_ge.sjdb_gtf_file == "-" {
            return Ok(Self {
                gtf_yes: false,
                ..Self::default()
            });
        }

        let dir_out = dir_out.as_ref();
        std::fs::create_dir_all(dir_out)?;

        if genome.chr_name_index.is_empty() {
            for (idx, chr) in genome.chr_name.iter().enumerate() {
                genome.chr_name_index.insert(chr.clone(), idx as u64);
            }
        }

        let reader =
            open_maybe_compressed(Path::new(&p.p_ge.sjdb_gtf_file)).with_context(|| {
                format!(
                    "FATAL error, could not open file pGe.sjdbGTFfile={}",
                    p.p_ge.sjdb_gtf_file
                )
            })?;

        let mut out = Self {
            gtf_yes: true,
            ..Self::default()
        };

        let mut transcript_id_number: HashMap<String, u64> = HashMap::new();
        let mut gene_id_number: HashMap<String, u64> = HashMap::new();
        let mut exon_lines_seen = 0usize;

        for line in reader.lines() {
            let line = line?;
            let Some(rec) = parse_exon_record(&line, p) else {
                continue;
            };
            exon_lines_seen += 1;

            let mut chr = rec.chr;
            if p.p_ge.sjdb_gtf_chr_prefix != "-" {
                chr = format!("{}{}", p.p_ge.sjdb_gtf_chr_prefix, chr);
            }

            let Some(&chr_idx_u64) = genome.chr_name_index.get(&chr) else {
                eprintln!(
                    "WARNING: while processing sjdbGTFfile={}: chromosome '{}' not found in Genome fasta files for line:\n{}",
                    p.p_ge.sjdb_gtf_file, chr, line
                );
                continue;
            };
            let chr_idx = chr_idx_u64 as usize;
            if rec.end > genome.chr_length[chr_idx] {
                eprintln!(
                    "WARNING: while processing sjdbGTFfile={}, line:\n{}\n exon end = {} is larger than the chromosome {} length = {} , will skip this exon",
                    p.p_ge.sjdb_gtf_file,
                    line,
                    rec.end,
                    chr,
                    genome.chr_length[chr_idx]
                );
                continue;
            }

            let tr_names = vec![p.p_ge.sjdb_gtf_tag_exon_parent_transcript.clone()];
            let ge_names = vec![p.p_ge.sjdb_gtf_tag_exon_parent_gene.clone()];
            let mut ex_attr = [
                extract_attr(&rec.attrs, &tr_names).unwrap_or_default(),
                extract_attr(&rec.attrs, &ge_names).unwrap_or_default(),
                extract_attr(&rec.attrs, &p.p_ge.sjdb_gtf_tag_exon_parent_gene_name)
                    .unwrap_or_default(),
                extract_attr(&rec.attrs, &p.p_ge.sjdb_gtf_tag_exon_parent_gene_type)
                    .unwrap_or_default(),
            ];

            if ex_attr[0].is_empty() {
                eprintln!(
                    "WARNING: while processing pGe.sjdbGTFfile={}: no transcript_id for line:\n{}",
                    p.p_ge.sjdb_gtf_file, line
                );
                ex_attr[0] = format!(
                    "tr_{}_{}_{}_{}",
                    chr,
                    rec.start,
                    rec.end,
                    out.exon_loci.len()
                );
            }
            if ex_attr[1].is_empty() {
                eprintln!(
                    "WARNING: while processing pGe.sjdbGTFfile={}: no gene_id for line:\n{}",
                    p.p_ge.sjdb_gtf_file, line
                );
                ex_attr[1] = "MissingGeneID".to_string();
            }
            if ex_attr[2].is_empty() {
                ex_attr[2] = ex_attr[1].clone();
            }
            if ex_attr[3].is_empty() {
                ex_attr[3] = "MissingGeneType".to_string();
            }

            let tr_idx = if let Some(&idx) = transcript_id_number.get(&ex_attr[0]) {
                idx
            } else {
                let idx = out.transcript_id.len() as u64;
                transcript_id_number.insert(ex_attr[0].clone(), idx);
                out.transcript_id.push(ex_attr[0].clone());
                out.transcript_strand.push(match rec.strand {
                    b'+' => 1,
                    b'-' => 2,
                    _ => 0,
                });
                idx
            };

            let ge_idx = if let Some(&idx) = gene_id_number.get(&ex_attr[1]) {
                idx
            } else {
                let idx = out.gene_id.len() as u64;
                gene_id_number.insert(ex_attr[1].clone(), idx);
                out.gene_id.push(ex_attr[1].clone());
                out.gene_attr.push([ex_attr[2].clone(), ex_attr[3].clone()]);
                idx
            };

            out.exon_loci.push([
                tr_idx,
                rec.start + genome.chr_start[chr_idx] - 1,
                rec.end + genome.chr_start[chr_idx] - 1,
                ge_idx,
            ]);
        }

        if exon_lines_seen == 0 {
            anyhow::bail!(
                "Fatal INPUT FILE error, no \"exon\" lines in the GTF file: {}\n\
                 Solution: check the formatting of the GTF file, it must contain some lines with \"exon\" in the 3rd column.\n\
                 Make sure the GTF file is unzipped.\n\
                 If exons are marked with a different word, use --sjdbGTFfeatureExon .",
                p.p_ge.sjdb_gtf_file
            );
        }
        if out.exon_loci.is_empty() {
            anyhow::bail!(
                "Fatal INPUT FILE error, no valid \"exon\" lines in the GTF file: {}\n\
                 Solution: check the formatting of the GTF file. One likely cause is the difference in chromosome naming between GTF and FASTA file.",
                p.p_ge.sjdb_gtf_file
            );
        }

        Ok(out)
    }

    pub fn transcript_gene_sj(
        &mut self,
        genome: &Genome,
        p: &Parameters,
        dir_out: impl AsRef<Path>,
        sjdb_loci: &mut SjdbLoci,
    ) -> Result<u64> {
        if !self.gtf_yes {
            return Ok(0);
        }

        self.exon_loci
            .sort_by(|a, b| a[EX_T].cmp(&b[EX_T]).then_with(|| a[EX_S].cmp(&b[EX_S])));
        let exon_n = self.exon_loci.len();
        let dir_out = dir_out.as_ref();
        std::fs::create_dir_all(dir_out)?;

        let mut exge_loci = vec![[0u64; EXGE_LEN]; exon_n];
        for (iex, exon) in self.exon_loci.iter().enumerate() {
            exge_loci[iex][EXGE_EX_START] = exon[EX_S];
            exge_loci[iex][EXGE_EX_END] = exon[EX_E];
            exge_loci[iex][EXGE_EX_STRAND] = self.transcript_strand[exon[EX_T] as usize] as u64;
            exge_loci[iex][EXGE_GE_ID] = exon[EX_G];
            exge_loci[iex][EXGE_TR_ID] = exon[EX_T];
        }
        exge_loci.sort();

        let mut exge_out = create_tab(dir_out.join("exonGeTrInfo.tab"))?;
        writeln!(exge_out, "{exon_n}")?;
        for row in &exge_loci {
            writeln!(
                exge_out,
                "{}\t{}\t{}\t{}\t{}",
                row[EXGE_EX_START],
                row[EXGE_EX_END],
                row[EXGE_EX_STRAND],
                row[EXGE_GE_ID],
                row[EXGE_TR_ID]
            )?;
        }
        exge_out.flush()?;

        let mut ge_out = create_tab(dir_out.join("geneInfo.tab"))?;
        writeln!(ge_out, "{}", self.gene_id.len())?;
        for (idx, gene) in self.gene_id.iter().enumerate() {
            writeln!(
                ge_out,
                "{}\t{}\t{}",
                gene, self.gene_attr[idx][0], self.gene_attr[idx][1]
            )?;
        }
        ge_out.flush()?;

        let mut extr_loci = vec![[0u64; EXTR_LEN]; exon_n];
        let mut trex1 = 0usize;
        for iex in 0..=exon_n {
            if iex == exon_n || self.exon_loci[iex][EX_T] != self.exon_loci[trex1][EX_T] {
                for iex1 in trex1..iex {
                    extr_loci[iex1][EXTR_TR_END] = self.exon_loci[iex - 1][EX_E];
                }
                if iex == exon_n {
                    break;
                }
                trex1 = iex;
            }
            extr_loci[iex][EXTR_TR_START] = self.exon_loci[trex1][EX_S];
            extr_loci[iex][EXTR_TR_ID] = self.exon_loci[iex][EX_T];
            extr_loci[iex][EXTR_EX_START] = self.exon_loci[iex][EX_S];
            extr_loci[iex][EXTR_EX_END] = self.exon_loci[iex][EX_E];
            extr_loci[iex][EXTR_GE_ID] = self.exon_loci[iex][EX_G];
        }
        extr_loci.sort_by(|a, b| a[..5].cmp(&b[..5]));

        let mut tr_out = create_tab(dir_out.join("transcriptInfo.tab"))?;
        writeln!(tr_out, "{}", self.transcript_id.len())?;
        let mut ex_out = create_tab(dir_out.join("exonInfo.tab"))?;
        writeln!(ex_out, "{exon_n}")?;

        let mut trid = extr_loci[0][EXTR_TR_ID];
        let mut trex = 0usize;
        let mut trstart = extr_loci[0][EXTR_TR_START];
        let mut trend = extr_loci[0][EXTR_TR_END];
        let mut exlen = 0u64;
        for iex in 0..=exon_n {
            if iex == exon_n || extr_loci[iex][EXTR_TR_ID] != trid {
                writeln!(
                    tr_out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    self.transcript_id[trid as usize],
                    extr_loci[iex - 1][EXTR_TR_START],
                    extr_loci[iex - 1][EXTR_TR_END],
                    trend,
                    self.transcript_strand[trid as usize] as u64,
                    iex - trex,
                    trex,
                    extr_loci[iex - 1][EXTR_GE_ID]
                )?;
                if iex == exon_n {
                    break;
                }
                trid = extr_loci[iex][EXTR_TR_ID];
                trstart = extr_loci[iex][EXTR_TR_START];
                trex = iex;
                trend = trend.max(extr_loci[iex - 1][EXTR_TR_END]);
                exlen = 0;
            }
            writeln!(
                ex_out,
                "{}\t{}\t{}",
                extr_loci[iex][EXTR_EX_START] - trstart,
                extr_loci[iex][EXTR_EX_END] - trstart,
                exlen
            )?;
            exlen += extr_loci[iex][EXTR_EX_END] - extr_loci[iex][EXTR_EX_START] + 1;
        }
        tr_out.flush()?;
        ex_out.flush()?;

        let mut sj_loci: Vec<[u64; 4]> = Vec::with_capacity(exon_n);
        let mut tr_idn = self.exon_loci[0][EX_T];
        for iex in 1..exon_n {
            if tr_idn == self.exon_loci[iex][EX_T] {
                let chr1 = genome.chr_bin
                    [(self.exon_loci[iex][EX_S] >> p.p_ge.g_chr_bin_nbits) as usize]
                    as usize;
                if self.exon_loci[iex][EX_S] <= self.exon_loci[iex - 1][EX_E] + 1 {
                    if self.exon_loci[iex][EX_S] <= self.exon_loci[iex - 1][EX_E] {
                        eprintln!(
                            "WARNING: while processing pGe.sjdbGTFfile={}: overlapping exons:\n{}\t{}\t{}\n{}\t{}\t{}",
                            p.p_ge.sjdb_gtf_file,
                            genome.chr_name[chr1],
                            self.exon_loci[iex - 1][EX_S] + 1 - genome.chr_start[chr1],
                            self.exon_loci[iex - 1][EX_E] + 1 - genome.chr_start[chr1],
                            genome.chr_name[chr1],
                            self.exon_loci[iex][EX_S] + 1 - genome.chr_start[chr1],
                            self.exon_loci[iex][EX_E] + 1 - genome.chr_start[chr1]
                        );
                    }
                } else {
                    sj_loci.push([
                        self.exon_loci[iex - 1][EX_E] + 1,
                        self.exon_loci[iex][EX_S] - 1,
                        self.transcript_strand[tr_idn as usize] as u64,
                        self.exon_loci[iex][EX_G] + 1,
                    ]);
                }
            } else {
                tr_idn = self.exon_loci[iex][EX_T];
            }
        }
        sj_loci.sort_by(|a, b| a[0].cmp(&b[0]).then_with(|| a[1].cmp(&b[1])));

        let strand_char = [b'.', b'+', b'-'];
        let sjdb_n1 = sjdb_loci.chr.len();
        sjdb_loci.gene.resize(sjdb_n1, BTreeSet::new());
        for (ii, row) in sj_loci.iter().enumerate() {
            let same_as_prev = ii > 0
                && row[0] == sj_loci[ii - 1][0]
                && row[1] == sj_loci[ii - 1][1]
                && row[2] == sj_loci[ii - 1][2];
            if !same_as_prev {
                let chr1 = genome.chr_bin[(row[0] >> p.p_ge.g_chr_bin_nbits) as usize] as usize;
                sjdb_loci.chr.push(genome.chr_name[chr1].clone());
                sjdb_loci.start.push(row[0] + 1 - genome.chr_start[chr1]);
                sjdb_loci.end.push(row[1] + 1 - genome.chr_start[chr1]);
                sjdb_loci.str_.push(strand_char[row[2] as usize]);
                sjdb_loci.gene.push(BTreeSet::from([row[3]]));
            } else {
                sjdb_loci.gene.last_mut().unwrap().insert(row[3]);
            }
        }

        let mut sjdb_list = create_tab(dir_out.join("sjdbList.fromGTF.out.tab"))?;
        for ii in sjdb_n1..sjdb_loci.chr.len() {
            write!(
                sjdb_list,
                "{}\t{}\t{}\t{}",
                sjdb_loci.chr[ii],
                sjdb_loci.start[ii],
                sjdb_loci.end[ii],
                sjdb_loci.str_[ii] as char
            )?;
            let mut genes = sjdb_loci.gene[ii].iter();
            if let Some(gene) = genes.next() {
                write!(sjdb_list, "\t{gene}")?;
            }
            for gene in genes {
                write!(sjdb_list, ",{gene}")?;
            }
            writeln!(sjdb_list)?;
        }
        sjdb_list.flush()?;

        sjdb_loci.priority.resize(sjdb_loci.chr.len(), 20);
        Ok((sjdb_loci.chr.len() - sjdb_n1) as u64)
    }
}

#[derive(Debug)]
struct ExonRecord {
    chr: String,
    start: u64,
    end: u64,
    strand: u8,
    attrs: String,
}

fn parse_exon_record(line: &str, p: &Parameters) -> Option<ExonRecord> {
    if line.starts_with('#') {
        return None;
    }
    let fields: Vec<&str> = line.splitn(9, '\t').collect();
    if fields.len() < 9 || fields[2] != p.p_ge.sjdb_gtf_feature_exon {
        return None;
    }
    Some(ExonRecord {
        chr: fields[0].to_string(),
        start: fields[3].parse().ok()?,
        end: fields[4].parse().ok()?,
        strand: fields[6].bytes().next().unwrap_or(b'.'),
        attrs: fields[8].to_string(),
    })
}

fn extract_attr(attrs: &str, names: &[String]) -> Option<String> {
    let normalized = normalize_attrs(attrs);
    let tokens: Vec<&str> = normalized.split_whitespace().collect();
    for name in names {
        for window in tokens.windows(2) {
            if window[0] == name {
                return Some(window[1].to_string());
            }
        }
    }
    None
}

fn normalize_attrs(attrs: &str) -> String {
    attrs
        .chars()
        .map(|ch| match ch {
            ';' | '=' | '\t' | '"' => ' ',
            other => other,
        })
        .collect()
}

fn create_tab(path: impl AsRef<Path>) -> Result<BufWriter<File>> {
    let path = path.as_ref();
    File::create(path)
        .map(BufWriter::new)
        .with_context(|| format!("could not create {}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;
    use std::process::Command;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn repo_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(3)
            .unwrap()
            .to_path_buf()
    }

    fn fixture_dir() -> PathBuf {
        repo_root().join("star-rs/tests/fixtures/gctest")
    }

    fn ref_bin() -> Option<PathBuf> {
        let candidate = repo_root().join("STAR/bin/Linux_x86_64/STAR");
        candidate.exists().then_some(candidate)
    }

    struct TempDir {
        path: PathBuf,
    }

    impl TempDir {
        fn new(prefix: &str) -> Self {
            let stamp = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let path = std::env::temp_dir().join(format!(
                "star-rs-{}-{}-{}",
                prefix,
                std::process::id(),
                stamp
            ));
            std::fs::create_dir_all(&path).unwrap();
            Self { path }
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.path);
        }
    }

    fn genome_generate_params(genome_dir: &Path, gtf: &Path) -> Parameters {
        let fixture = fixture_dir();
        let args = vec![
            "STAR".to_string(),
            "--runMode".to_string(),
            "genomeGenerate".to_string(),
            "--runThreadN".to_string(),
            "1".to_string(),
            "--genomeDir".to_string(),
            genome_dir.display().to_string(),
            "--genomeFastaFiles".to_string(),
            fixture.join("chr.fa").display().to_string(),
            "--genomeSAindexNbases".to_string(),
            "5".to_string(),
            "--sjdbGTFfile".to_string(),
            gtf.display().to_string(),
            "--sjdbOverhang".to_string(),
            "49".to_string(),
            "--outFileNamePrefix".to_string(),
            format!("{}/", genome_dir.display()),
        ];
        let mut p = Parameters::new();
        p.parse_cli(&args).unwrap();
        p.command_line_full = args.join(" ");
        p
    }

    fn run_ref_genome_generate(genome_dir: &Path, gtf: &Path) {
        let Some(bin) = ref_bin() else {
            return;
        };
        let fixture = fixture_dir();
        let output = Command::new(bin)
            .args([
                "--runMode",
                "genomeGenerate",
                "--runThreadN",
                "1",
                "--genomeDir",
                &genome_dir.display().to_string(),
                "--genomeFastaFiles",
                &fixture.join("chr.fa").display().to_string(),
                "--genomeSAindexNbases",
                "5",
                "--sjdbGTFfile",
                &gtf.display().to_string(),
                "--sjdbOverhang",
                "49",
                "--outFileNamePrefix",
                &format!("{}/", genome_dir.display()),
            ])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "reference genomeGenerate failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    fn read_file(path: &Path) -> String {
        std::fs::read_to_string(path).unwrap()
    }

    #[test]
    fn gtf_new_parses_custom_attribute_names() {
        let dir = TempDir::new("gtf-attrs");
        let fasta = dir.path.join("tiny.fa");
        let gtf = dir.path.join("tiny.gtf");
        std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
        std::fs::write(
            &gtf,
            "chr1\tmk\texon\t2\t6\t.\t+\t.\tParentTx=tx1; ParentGene=gene1; Name=Gene1; Kind=protein_coding;\n",
        )
        .unwrap();

        let args = vec![
            "STAR".to_string(),
            "--runMode".to_string(),
            "genomeGenerate".to_string(),
            "--genomeDir".to_string(),
            dir.path.join("Genome").display().to_string(),
            "--genomeFastaFiles".to_string(),
            fasta.display().to_string(),
            "--genomeSAindexNbases".to_string(),
            "3".to_string(),
            "--sjdbGTFfile".to_string(),
            gtf.display().to_string(),
            "--sjdbOverhang".to_string(),
            "5".to_string(),
            "--sjdbGTFtagExonParentTranscript".to_string(),
            "ParentTx".to_string(),
            "--sjdbGTFtagExonParentGene".to_string(),
            "ParentGene".to_string(),
            "--sjdbGTFtagExonParentGeneName".to_string(),
            "Name".to_string(),
            "--sjdbGTFtagExonParentGeneType".to_string(),
            "Kind".to_string(),
        ];
        let mut p = Parameters::new();
        p.parse_cli(&args).unwrap();

        let mut genome = star_genome::genome_generate(&mut p).unwrap();
        genome.chr_bin_fill();
        let mut sjdb_loci = SjdbLoci::new();
        let gtf = Gtf::new(&mut genome, &p, &dir.path, &mut sjdb_loci).unwrap();
        assert!(gtf.gtf_yes);
        assert_eq!(gtf.transcript_id, vec!["tx1"]);
        assert_eq!(gtf.gene_id, vec!["gene1"]);
        assert_eq!(
            gtf.gene_attr,
            vec![["Gene1".to_string(), "protein_coding".to_string()]]
        );
        assert_eq!(gtf.exon_loci.len(), 1);
    }

    #[test]
    fn transcript_gene_sj_matches_reference_gctest_files() {
        let Some(_) = ref_bin() else {
            eprintln!("STAR reference binary missing, skipping");
            return;
        };

        let rs_dir = TempDir::new("gtf-rs");
        let ref_dir = TempDir::new("gtf-ref");
        let fixture = fixture_dir();
        let gtf_path = fixture.join("ann.gtf");

        let mut p = genome_generate_params(&rs_dir.path, &gtf_path);
        let mut genome = star_genome::genome_generate(&mut p).unwrap();
        genome.chr_bin_fill();
        let mut sjdb_loci = SjdbLoci::new();
        let mut gtf = Gtf::new(&mut genome, &p, &rs_dir.path, &mut sjdb_loci).unwrap();
        gtf.transcript_gene_sj(&genome, &p, &rs_dir.path, &mut sjdb_loci)
            .unwrap();

        run_ref_genome_generate(&ref_dir.path, &gtf_path);

        for name in [
            "geneInfo.tab",
            "transcriptInfo.tab",
            "exonInfo.tab",
            "exonGeTrInfo.tab",
            "sjdbList.fromGTF.out.tab",
        ] {
            assert_eq!(
                read_file(&rs_dir.path.join(name)),
                read_file(&ref_dir.path.join(name)),
                "mismatch in {name}"
            );
        }
        assert!(sjdb_loci.priority.iter().all(|&priority| priority == 20));
    }
}
