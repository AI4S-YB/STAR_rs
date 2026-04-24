use std::io::Write;
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

fn run_cmd(args: &[String]) {
    let bin = rust_bin();
    let out = Command::new(&bin).args(args).output().unwrap();
    assert!(
        out.status.success(),
        "command failed: {}\nstdout:\n{}\nstderr:\n{}",
        bin.display(),
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}

fn write_gzip(src: &Path, dst: &Path) {
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    enc.write_all(&std::fs::read(src).unwrap()).unwrap();
    std::fs::write(dst, enc.finish().unwrap()).unwrap();
}

fn write_xz(src: &Path, dst: &Path) {
    let mut enc = xz2::write::XzEncoder::new(Vec::new(), 6);
    enc.write_all(&std::fs::read(src).unwrap()).unwrap();
    std::fs::write(dst, enc.finish().unwrap()).unwrap();
}

fn write_bzip2(src: &Path, dst: &Path) {
    let mut enc = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    enc.write_all(&std::fs::read(src).unwrap()).unwrap();
    std::fs::write(dst, enc.finish().unwrap()).unwrap();
}

fn genome_generate_args(genome_dir: &Path, fasta: &Path, gtf: &Path) -> Vec<String> {
    vec![
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--runThreadN".to_string(),
        "1".to_string(),
        "--genomeDir".to_string(),
        genome_dir.display().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.display().to_string(),
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

fn align_args(genome_dir: &Path, read: &Path, out_prefix: &Path) -> Vec<String> {
    vec![
        "--runMode".to_string(),
        "alignReads".to_string(),
        "--runThreadN".to_string(),
        "1".to_string(),
        "--genomeDir".to_string(),
        genome_dir.display().to_string(),
        "--readFilesIn".to_string(),
        read.display().to_string(),
        "--quantMode".to_string(),
        "GeneCounts".to_string(),
        "--outFileNamePrefix".to_string(),
        format!("{}/", out_prefix.display()),
    ]
}

fn sam_body(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('@'))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn compressed_genome_gtf_and_reads_match_plain_inputs() {
    let fixture = fixture_dir();
    let input_dir = TempDir::new("compressed-inputs");
    let plain_genome = TempDir::new("plain-genome");
    let compressed_genome = TempDir::new("compressed-genome");
    let plain_align = TempDir::new("plain-align");
    let compressed_align = TempDir::new("compressed-align");

    let fasta_gz = input_dir.path.join("chr.fa.gz");
    let gtf_xz = input_dir.path.join("ann.gtf.xz");
    let reads_bz = input_dir.path.join("r.fq.bz");
    write_gzip(&fixture.join("chr.fa"), &fasta_gz);
    write_xz(&fixture.join("ann.gtf"), &gtf_xz);
    write_bzip2(&fixture.join("r.fq"), &reads_bz);

    run_cmd(&genome_generate_args(
        &plain_genome.path,
        &fixture.join("chr.fa"),
        &fixture.join("ann.gtf"),
    ));
    run_cmd(&genome_generate_args(
        &compressed_genome.path,
        &fasta_gz,
        &gtf_xz,
    ));

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
            std::fs::read(compressed_genome.path.join(name)).unwrap(),
            std::fs::read(plain_genome.path.join(name)).unwrap(),
            "mismatch in {name}"
        );
    }

    run_cmd(&align_args(
        &compressed_genome.path,
        &fixture.join("r.fq"),
        &plain_align.path,
    ));
    run_cmd(&align_args(
        &compressed_genome.path,
        &reads_bz,
        &compressed_align.path,
    ));

    assert_eq!(
        sam_body(&compressed_align.path.join("Aligned.out.sam")),
        sam_body(&plain_align.path.join("Aligned.out.sam"))
    );
    assert_eq!(
        std::fs::read(compressed_align.path.join("SJ.out.tab")).unwrap(),
        std::fs::read(plain_align.path.join("SJ.out.tab")).unwrap()
    );
    assert_eq!(
        std::fs::read(compressed_align.path.join("ReadsPerGene.out.tab")).unwrap(),
        std::fs::read(plain_align.path.join("ReadsPerGene.out.tab")).unwrap()
    );
}
