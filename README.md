# star-rs

`STAR` RNA-seq aligner — 1:1 Rust port of
[`alexdobin/STAR`](https://github.com/alexdobin/STAR) (C++ reference lives at
`../STAR/`).

Goal: **byte-exact** / **bit-exact** parity with the C++ binary for the
supported subset, verified by the end-to-end regression script in
[`tests/e2e.sh`](tests/e2e.sh).

## Status

| Milestone | Scope | State |
|-----------|-------|-------|
| M1 | Workspace scaffolding, crate graph, CLI skeleton | done |
| M2 | `--runMode genomeGenerate` (bit-exact `Genome`/`SA`/`SAindex`/`chr*`) | done |
| M3 | `--runMode alignReads` SE (single thread) | done |
| M4 | PE + `--runThreadN>1` with deterministic chunked output | done |
| M5 | `--sjdbFileChrStartEnd` + `--twopassMode Basic` | done |
| M6 | Chimeric detection (legacy Old path) byte-exact | done |
| M7.1–M7.3 | `--quantMode GeneCounts` → `ReadsPerGene.out.tab` byte-exact | done |
| M7.4 | `--outSAMtype BAM Unsorted` **minimum viable** (semantic parity, not byte-exact) | done |
| M7-parked | remaining BAM/wig/dedup features (see `docs/roadmap.md`) | pending |
| M8-tier1 | Internal regression on Mp/chr1 × 10k PE reads (byte-exact) | done |
| M8-tier2 | Full *M. polymorpha* genome × 10k PE reads | pending |
| M8-tier3 | Full *M. polymorpha* genome × MP-S1 library (gzipped reads) | blocked on `readFilesCommand` |

Regression coverage: 29 cases in `tests/e2e.sh` (all green) plus the
Mp/chr1 tier-1 matrix and checked-in fixture regressions described in
[`docs/regression-notes.md`](docs/regression-notes.md).

## Layout

```
star-rs/
├── Cargo.toml                # workspace
├── crates/
│   ├── star-core/            # shared types, constants, Genome version
│   ├── star-params/          # CLI parsing + `Parameters` finalize()
│   ├── star-genome/          # Genome / SA / SAindex / sjdb state
│   ├── star-io/              # SAM headers, SJ writer, fastq readers
│   ├── star-align/           # ReadAlign pipeline, chunk + thread orchestration
│   ├── star-chimeric/        # chimeric detection (Old path)
│   ├── star-quant/           # Transcriptome loader, GeneCounts
│   ├── star-bam/             # SAM→BAM converter (noodles, minimum viable)
│   ├── star-sjdb/            # sjdb prepare / build index / 2-pass glue
│   ├── star-stats/           # Log.final.out counters
│   ├── opal-sys/             # FFI to opal.o
│   └── star-cli/             # `star` binary (run_mapping_pass, runMode dispatch)
└── tests/
    ├── e2e.sh                # driver: run Rust + C++ STAR, diff outputs
    └── fixtures/             # tiny / sjtest / chimtest / gctest
```

## Building

```bash
cargo build --release -p star-cli
./target/release/star --help
```

## Testing

```bash
./tests/e2e.sh
# STAR_REF=/path/STAR ./tests/e2e.sh    override C++ reference
# STAR_RS=target/release/star ./tests/e2e.sh
# KEEP=1 ./tests/e2e.sh                 keep scratch dir
```

## Current validated coverage

- `genomeGenerate`:
  - `Genome` / `SA` / `SAindex` / `chr*.txt` byte-exact on current checked-in fixtures
  - `--sjdbGTFfile` supported, including
    `geneInfo.tab` / `transcriptInfo.tab` / `exonInfo.tab` /
    `exonGeTrInfo.tab` / `sjdbList.fromGTF.out.tab`
- `alignReads`:
  - SE + PE
  - deterministic multithreaded output
  - `--sjdbFileChrStartEnd`
  - `--twopassMode Basic`
  - chimeric Old path
  - `--quantMode GeneCounts`
  - `--outSAMtype BAM Unsorted` minimum viable output
- Byte-exact checks currently cover:
  - `Aligned.out.sam` body
  - `SJ.out.tab`
  - `ReadsPerGene.out.tab`
  - `Log.final.out` stats lines on the checked-in `tests/` fixture

## Performance

On the checked-in root `tests/` fixture (`chr1.fa + chr1.gtf`,
`reads_10k_{1,2}.fq`, `--runThreadN 4`, `--quantMode GeneCounts`):

- `alignReads`: Rust averaged `14.98 s` across 3 runs vs `30.44 s`
  for C++ STAR, with peak RSS `494.7 MiB` vs `745.0 MiB`.
- `genomeGenerate`: Rust is still slower/heavier on this fixture
  (`41.90 s`, `1128.9 MiB`) than C++ STAR (`13.10 s`, `920.6 MiB`).

Full measurement details are in [docs/regression-notes.md](docs/regression-notes.md).

## Documentation

- [docs/roadmap.md](docs/roadmap.md)
  M7-parked items, M8 plan, internal regression tiers, and remaining work.
- [docs/regression-notes.md](docs/regression-notes.md)
  Detailed debug history, M8-tier1 notes, Ath regression notes, and performance snapshots.
- [docs/star_rust_port_plan_3ec86fb7.plan.md](docs/star_rust_port_plan_3ec86fb7.plan.md)
  Original detailed port plan.

## License

MIT, matching the upstream C++ project.
