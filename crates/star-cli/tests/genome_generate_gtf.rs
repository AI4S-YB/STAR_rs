use std::path::{Path, PathBuf};
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

fn rust_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_star"))
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

fn run_cmd(bin: &Path, args: &[String]) {
    let out = Command::new(bin).args(args).output().unwrap();
    assert!(
        out.status.success(),
        "command failed: {}\nstdout:\n{}\nstderr:\n{}",
        bin.display(),
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}

fn genome_generate_args(genome_dir: &Path, gtf: &Path) -> Vec<String> {
    let fixture = fixture_dir();
    vec![
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
    ]
}

fn align_gene_counts_args(genome_dir: &Path, out_prefix: &Path) -> Vec<String> {
    let fixture = fixture_dir();
    vec![
        "--runMode".to_string(),
        "alignReads".to_string(),
        "--runThreadN".to_string(),
        "1".to_string(),
        "--genomeDir".to_string(),
        genome_dir.display().to_string(),
        "--readFilesIn".to_string(),
        fixture.join("r.fq").display().to_string(),
        "--quantMode".to_string(),
        "GeneCounts".to_string(),
        "--outFileNamePrefix".to_string(),
        format!("{}/", out_prefix.display()),
    ]
}

fn normalized_genome_parameters(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn genome_generate_with_gtf_matches_reference_and_reloads_for_gene_counts() {
    let Some(ref_star) = ref_bin() else {
        eprintln!("STAR reference binary missing, skipping");
        return;
    };

    let fixture = fixture_dir();
    let gtf = fixture.join("ann.gtf");

    let rs_genome = TempDir::new("rs-genome");
    let ref_genome = TempDir::new("ref-genome");
    let rs_align = TempDir::new("rs-align");
    let ref_align = TempDir::new("ref-align");

    run_cmd(&rust_bin(), &genome_generate_args(&rs_genome.path, &gtf));
    run_cmd(&ref_star, &genome_generate_args(&ref_genome.path, &gtf));

    for name in [
        "Genome",
        "SA",
        "SAindex",
        "chrName.txt",
        "chrStart.txt",
        "chrLength.txt",
        "chrNameLength.txt",
        "sjdbInfo.txt",
        "sjdbList.out.tab",
        "geneInfo.tab",
        "transcriptInfo.tab",
        "exonInfo.tab",
        "exonGeTrInfo.tab",
        "sjdbList.fromGTF.out.tab",
    ] {
        assert_eq!(
            std::fs::read(rs_genome.path.join(name)).unwrap(),
            std::fs::read(ref_genome.path.join(name)).unwrap(),
            "mismatch in {name}"
        );
    }

    assert_eq!(
        normalized_genome_parameters(&rs_genome.path.join("genomeParameters.txt")),
        normalized_genome_parameters(&ref_genome.path.join("genomeParameters.txt"))
    );

    run_cmd(
        &rust_bin(),
        &align_gene_counts_args(&rs_genome.path, &rs_align.path),
    );
    run_cmd(
        &ref_star,
        &align_gene_counts_args(&ref_genome.path, &ref_align.path),
    );

    assert_eq!(
        std::fs::read(rs_align.path.join("ReadsPerGene.out.tab")).unwrap(),
        std::fs::read(ref_align.path.join("ReadsPerGene.out.tab")).unwrap()
    );
}
