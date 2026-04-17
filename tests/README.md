# star-rs regression tests

Byte-exact end-to-end tests comparing the Rust STAR port against the C++
reference binary.

## Layout

```
tests/
├── e2e.sh                 # driver: run every case, cmp vs C++ STAR
├── fixtures/
│   ├── tiny/              # 4 KB chr1 + 10 reads (SE + PE)
│   │   ├── chr1.fa
│   │   ├── reads.fq
│   │   ├── reads_r1.fq
│   │   └── reads_r2.fq
│   └── sjtest/            # 1.2 KB ref with canonical GT-AG intron
│       ├── chr.fa
│       ├── r.fq
│       ├── sj.tab         # --sjdbFileChrStartEnd input
│       └── mk_sj.py       # deterministic regeneration script
```

## Running

```bash
cargo build --release -p star-cli
./tests/e2e.sh
```

Environment knobs:

- `STAR_REF=/path/to/STAR` — override C++ STAR binary (default: `../STAR/bin/Linux_x86_64/STAR`)
- `STAR_RS=target/release/star` — override Rust binary
- `KEEP=1` — keep scratch directory for inspection

## What's covered

| # | Case | Checked outputs |
|---|------|-----------------|
| 1 | `genomeGenerate` on tiny/chr1.fa | `Genome`, `SA`, `SAindex`, `chr*.txt` bit-identical |
| 2 | `alignReads` SE, `runThreadN=1` | `Aligned.out.sam` body, `SJ.out.tab` |
| 3 | `alignReads` PE, `runThreadN=4` | `Aligned.out.sam` body, `SJ.out.tab` |
| 4 | `alignReads --sjdbFileChrStartEnd` | genome bit-exact; `sjdbInfo.txt`, `sjdbList.out.tab` |
| 5 | `alignReads --twopassMode Basic` | SAM body, final `SJ.out.tab`, pass1 `SJ.out.tab`, `sjdbInfo.txt`, `sjdbList.out.tab` |

Only SAM bodies are compared (not `@HD`/`@PG` headers) because those embed
the CLI invocation and version string which differ between binaries.

## Adding a new case

1. Drop the fixture under `tests/fixtures/<name>/`.
2. Append a `[N] …` block to `e2e.sh` following the existing pattern
   (run both binaries into separate dirs, `cmp` selected files via
   `cmp_or` or `cmp_sam_body`).
