# Compressed Text Inputs

**Date:** 2026-04-24
**Status:** implemented in v0.3.1

## Goal

STAR-rs should accept common compressed text files directly for the text
input surfaces used by genome generation, annotation loading, splice-junction
loading, and FASTQ/FASTA read alignment. Users should not need to pre-stage
large `.fq.gz` files through `zcat` for the supported formats.

## Supported Formats

Compression is inferred from the file name suffix:

| Suffix | Decoder |
| ------ | ------- |
| `.gz`, `.gzip` | gzip |
| `.xz` | xz |
| `.bz`, `.bz2`, `.bzip2` | bzip2 |
| anything else | plain text |

Detection is extension-based by design. Plain files with non-compression
suffixes remain byte-for-byte plain text inputs.

## Input Surfaces

The shared opener is used for:

- `--readFilesIn`
- `--genomeFastaFiles`
- `--sjdbGTFfile`
- `--sjdbFileChrStartEnd`

This covers FASTQ, FASTA, GTF/GFF-like tabular annotation, and STAR splice
junction text inputs.

## Non-Goals

- `--readFilesCommand` shell-pipe compatibility is not implemented by this
  feature. It remains useful for arbitrary commands such as `samtools view`.
- BAM/SAM inputs that require external conversion are outside this feature.
- Compression is not auto-detected by magic bytes.

## Validation

- `cargo test -p star-core compression::tests`
- `cargo test -p star-cli --test compressed_inputs`

The CLI regression creates compressed FASTA, GTF, and FASTQ fixtures and
compares the generated genome sidecar files plus alignment outputs against
the same plain-text inputs.
