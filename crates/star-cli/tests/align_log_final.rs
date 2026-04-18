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

fn root_test_dir() -> PathBuf {
    repo_root().join("tests")
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

fn cpp_genome_generate_args(genome_dir: &Path) -> Vec<String> {
    let fixture = root_test_dir();
    vec![
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--runThreadN".to_string(),
        "4".to_string(),
        "--genomeDir".to_string(),
        genome_dir.display().to_string(),
        "--genomeFastaFiles".to_string(),
        fixture.join("chr1.fa").display().to_string(),
        "--genomeSAindexNbases".to_string(),
        "11".to_string(),
        "--sjdbGTFfile".to_string(),
        fixture.join("chr1.gtf").display().to_string(),
        "--sjdbOverhang".to_string(),
        "149".to_string(),
        "--outFileNamePrefix".to_string(),
        format!("{}/", genome_dir.display()),
    ]
}

fn align_args(genome_dir: &Path, out_prefix: &Path) -> Vec<String> {
    let fixture = root_test_dir();
    vec![
        "--runMode".to_string(),
        "alignReads".to_string(),
        "--runThreadN".to_string(),
        "4".to_string(),
        "--genomeDir".to_string(),
        genome_dir.display().to_string(),
        "--readFilesIn".to_string(),
        fixture.join("reads_10k_1.fq").display().to_string(),
        fixture.join("reads_10k_2.fq").display().to_string(),
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

fn normalized_log_final(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(4)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn align_log_final_matches_reference_stats_on_cpp_index() {
    let Some(ref_star) = ref_bin() else {
        eprintln!("STAR reference binary missing, skipping");
        return;
    };

    let genome_dir = TempDir::new("cpp-index");
    let cpp_align = TempDir::new("cpp-align");
    let rs_align = TempDir::new("rs-align");

    run_cmd(&ref_star, &cpp_genome_generate_args(&genome_dir.path));
    run_cmd(&ref_star, &align_args(&genome_dir.path, &cpp_align.path));
    run_cmd(&rust_bin(), &align_args(&genome_dir.path, &rs_align.path));

    assert_eq!(
        std::fs::read(rs_align.path.join("ReadsPerGene.out.tab")).unwrap(),
        std::fs::read(cpp_align.path.join("ReadsPerGene.out.tab")).unwrap()
    );
    assert_eq!(
        std::fs::read(rs_align.path.join("SJ.out.tab")).unwrap(),
        std::fs::read(cpp_align.path.join("SJ.out.tab")).unwrap()
    );
    assert_eq!(
        sam_body(&rs_align.path.join("Aligned.out.sam")),
        sam_body(&cpp_align.path.join("Aligned.out.sam"))
    );
    assert_eq!(
        normalized_log_final(&rs_align.path.join("Log.final.out")),
        normalized_log_final(&cpp_align.path.join("Log.final.out"))
    );
}
