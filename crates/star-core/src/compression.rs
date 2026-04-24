//! Shared text-input opening for plain and compressed STAR inputs.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Compression inferred from the input path suffix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionFormat {
    Plain,
    Gzip,
    Xz,
    Bzip2,
}

/// Infer compression from the file name. The check is intentionally
/// extension-based so plain FASTA/FASTQ/GTF parsing still sees the original
/// bytes when users choose non-compression suffixes.
pub fn compression_format_from_path(path: impl AsRef<Path>) -> CompressionFormat {
    let name = path
        .as_ref()
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    if name.ends_with(".gz") || name.ends_with(".gzip") {
        CompressionFormat::Gzip
    } else if name.ends_with(".xz") {
        CompressionFormat::Xz
    } else if name.ends_with(".bz") || name.ends_with(".bz2") || name.ends_with(".bzip2") {
        CompressionFormat::Bzip2
    } else {
        CompressionFormat::Plain
    }
}

/// Open a file as buffered text input, transparently decoding supported
/// compression formats based on the path suffix.
pub fn open_maybe_compressed(path: impl AsRef<Path>) -> io::Result<Box<dyn BufRead + Send>> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let reader: Box<dyn BufRead + Send> = match compression_format_from_path(path) {
        CompressionFormat::Plain => Box::new(BufReader::new(file)),
        CompressionFormat::Gzip => {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
        }
        CompressionFormat::Xz => Box::new(BufReader::new(xz2::read::XzDecoder::new(file))),
        CompressionFormat::Bzip2 => Box::new(BufReader::new(bzip2::read::BzDecoder::new(file))),
    };
    Ok(reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Write};
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    struct TempDir {
        path: PathBuf,
    }

    impl TempDir {
        fn new() -> Self {
            let stamp = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let path = std::env::temp_dir().join(format!(
                "star-rs-compression-{}-{}",
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

    fn read_to_string(path: &Path) -> io::Result<String> {
        let mut reader = open_maybe_compressed(path)?;
        let mut out = String::new();
        reader.read_to_string(&mut out)?;
        Ok(out)
    }

    #[test]
    fn detects_supported_suffixes() {
        assert_eq!(
            compression_format_from_path("reads.fq.gz"),
            CompressionFormat::Gzip
        );
        assert_eq!(
            compression_format_from_path("genome.FA.XZ"),
            CompressionFormat::Xz
        );
        assert_eq!(
            compression_format_from_path("ann.gtf.bz2"),
            CompressionFormat::Bzip2
        );
        assert_eq!(
            compression_format_from_path("ann.gtf.bz"),
            CompressionFormat::Bzip2
        );
        assert_eq!(
            compression_format_from_path("reads.fq"),
            CompressionFormat::Plain
        );
    }

    #[test]
    fn opens_plain_gzip_xz_and_bzip2() -> io::Result<()> {
        let dir = TempDir::new();
        let text = "alpha\nbeta\n";

        let plain = dir.path.join("plain.txt");
        std::fs::write(&plain, text)?;
        assert_eq!(read_to_string(&plain)?, text);

        let gzip = dir.path.join("reads.fq.gz");
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(text.as_bytes())?;
        std::fs::write(&gzip, enc.finish()?)?;
        assert_eq!(read_to_string(&gzip)?, text);

        let xz = dir.path.join("genome.fa.xz");
        let mut enc = xz2::write::XzEncoder::new(Vec::new(), 6);
        enc.write_all(text.as_bytes())?;
        std::fs::write(&xz, enc.finish()?)?;
        assert_eq!(read_to_string(&xz)?, text);

        let bzip2 = dir.path.join("ann.gtf.bz2");
        let mut enc = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
        enc.write_all(text.as_bytes())?;
        std::fs::write(&bzip2, enc.finish()?)?;
        assert_eq!(read_to_string(&bzip2)?, text);

        Ok(())
    }
}
