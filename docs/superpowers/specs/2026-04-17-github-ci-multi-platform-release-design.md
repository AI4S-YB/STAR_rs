# GitHub CI & Multi-Platform Release — Design

**Date:** 2026-04-17
**Repository:** https://github.com/AI4S-YB/STAR_rs
**Status:** Approved (pending user review of this document)

## Goal

Add GitHub Actions workflows that:

1. Run lint/format/test/MSRV checks on every PR and push to `main`.
2. Build and publish multi-platform release binaries when a `v*` tag is pushed (or on manual dispatch).

## Non-Goals

- No publishing to crates.io.
- No Docker image builds.
- No benchmark tracking.
- No code coverage reporting.

## Target Platforms (release binaries)

| Target triple | GitHub runner | Archive format |
|---|---|---|
| `x86_64-unknown-linux-gnu` | `ubuntu-latest` | `.tar.gz` |
| `x86_64-unknown-linux-musl` | `ubuntu-latest` | `.tar.gz` |
| `x86_64-apple-darwin` | `macos-13` | `.tar.gz` |
| `aarch64-apple-darwin` | `macos-14` | `.tar.gz` |
| `x86_64-pc-windows-msvc` | `windows-latest` | `.zip` |

Note on AVX2: `opal-sys/build.rs` compiles `opal.cpp` with `-mavx2` only when the upstream `STAR/source/opal` tree exists next to the workspace. CI won't have that tree, so `opal-sys` builds as a stub on every target. This is intentional and already handled by `build.rs`.

## Architecture

Two workflow files under `.github/workflows/`:

```
.github/
└── workflows/
    ├── ci.yml        # fmt / clippy / test / MSRV — on PR and push to main
    └── release.yml   # multi-platform build + publish — on v* tag or manual dispatch
```

Rationale for two files (not one, not three):

- **Not one** — CI and Release have different triggers, permissions (`contents: write` only for Release), and cadence. One file with nested `if` guards becomes unreadable.
- **Not three** — MSRV is cheap (`cargo check` only) and logically part of CI. Splitting it out multiplies checkout/cache overhead and hides the overall green/red state.

## CI workflow (`ci.yml`)

**Triggers:**
- `pull_request` targeting `main`
- `push` to `main`

**Concurrency:** `group: ci-${{ github.ref }}`, `cancel-in-progress: true` — a new push to the same PR/branch cancels the previous run.

**Jobs (parallel):**

| Job | Runner | Command |
|---|---|---|
| `fmt` | `ubuntu-latest` | `cargo fmt --all -- --check` |
| `clippy` | `ubuntu-latest` | `cargo clippy --workspace --all-targets -- -D warnings` |
| `test` | matrix: `ubuntu-latest` / `macos-latest` / `windows-latest` | `cargo test --workspace` |
| `msrv` | `ubuntu-latest`, rust `1.75` | `cargo check --workspace --all-targets` |

**Shared setup for each job:**
- `actions/checkout@v4`
- `dtolnay/rust-toolchain@stable` (or `@1.75` for `msrv`)
- `Swatinem/rust-cache@v2` with job-scoped key

**Matrix for `test`:**
```yaml
strategy:
  fail-fast: false
  matrix:
    os: [ubuntu-latest, macos-latest, windows-latest]
```

**MSRV rationale:** `Cargo.toml` declares `rust-version = "1.75"`. Running `cargo check` (not `test`) on 1.75 prevents silently raising the MSRV when someone uses a newer API, without requiring the test harness itself to stay on 1.75 (a common convention).

## Release workflow (`release.yml`)

**Triggers:**
- `push` on tags matching `v*`
- `workflow_dispatch` with input `tag` (default `v0.0.0-dev`) for manual/preview runs

**Permissions:** `contents: write` (needed to create/edit releases and upload assets).

**Jobs:**

### Job 1 — `create-release`

Runs once. Creates a **draft** GitHub Release for the tag (if not already present).

```
if gh release view "$TAG" >/dev/null 2>&1; then
  echo "Release $TAG already exists, skipping create"
else
  gh release create "$TAG" --draft --generate-notes --title "$TAG"
fi
```

`--generate-notes` auto-populates the body from commit messages since the previous tag. The `gh release view` check makes the step safely re-runnable (e.g., if a build matrix entry failed and the workflow is re-dispatched on the same tag).

### Job 2 — `build` (matrix)

Depends on `create-release`. `fail-fast: false` — one target failure does not cancel siblings.

Matrix:
```yaml
include:
  - target: x86_64-unknown-linux-gnu
    os: ubuntu-latest
  - target: x86_64-unknown-linux-musl
    os: ubuntu-latest
    apt: musl-tools
  - target: x86_64-apple-darwin
    os: macos-13
  - target: aarch64-apple-darwin
    os: macos-14
  - target: x86_64-pc-windows-msvc
    os: windows-latest
```

Each job:
1. `actions/checkout@v4`
2. `dtolnay/rust-toolchain@stable` with `targets: ${{ matrix.target }}`
3. `Swatinem/rust-cache@v2` keyed on target
4. (Linux-musl only) `sudo apt-get install -y ${{ matrix.apt }}`
5. `taiki-e/upload-rust-binary-action@v1` with:
   - `bin: star`
   - `target: ${{ matrix.target }}`
   - `archive: star-$tag-$target`
   - `checksum: sha256`
   - `tar: unix`
   - `zip: windows`

The action builds with `cargo build --release --bin star --target <triple>`, packages the artifact, generates `.sha256`, and uploads to the release associated with the current tag.

### Job 3 — `publish-release`

`needs: [create-release, build]`. Transitions the draft to a published release only if all build matrix jobs succeeded:

```
gh release edit "$TAG" --draft=false
```

If any `build` matrix entry fails, the release stays in draft state. The maintainer can investigate, fix, delete the draft, and re-run.

## Release asset naming

```
star-v0.1.0-x86_64-unknown-linux-gnu.tar.gz
star-v0.1.0-x86_64-unknown-linux-gnu.tar.gz.sha256
star-v0.1.0-x86_64-unknown-linux-musl.tar.gz
star-v0.1.0-x86_64-apple-darwin.tar.gz
star-v0.1.0-aarch64-apple-darwin.tar.gz
star-v0.1.0-x86_64-pc-windows-msvc.zip
```

Each archive contains a single `star` (or `star.exe`) binary at the root, plus `LICENSE` and `README.md` auto-included by `taiki-e/upload-rust-binary-action`.

## Secondary changes

Update `Cargo.toml`:
- `repository = "https://github.com/AI4S-YB/STAR_rs"` (currently placeholder `example/star-rs`)

No `README.md` badges in the initial change — can be added later once workflows are green.

## Security

- All third-party actions pinned to major version tags (`@v1`, `@v2`, `@v4`), matching GitHub's minimum recommended hygiene. Optional tightening to commit SHAs is deferred.
- No secrets other than auto-provided `GITHUB_TOKEN`.
- No caching of credentials across jobs.

## Testing the workflows

- **CI** — naturally exercised by the first PR that adds the workflow files.
- **Release** — two options:
  1. Push a throwaway tag like `v0.0.0-test1`, observe the run, delete the draft and the tag.
  2. Use `workflow_dispatch` with a custom `tag` input; the action will create a draft release on that tag.

Either way, verify:
- All five matrix jobs succeed.
- Six archive files + six `.sha256` files appear on the release.
- Downloading one archive, extracting, and running `./star --help` on the corresponding OS works.
- Release auto-transitions from draft to published.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Windows `cargo test` flaky due to path/newline differences | `fail-fast: false` in CI test matrix; triage and fix or narrow the target |
| A target's toolchain install changes silently | Pinned action major versions; rust-cache will rebuild on change |
| Clippy warnings from a new Rust release break main | Pinning clippy to a specific toolchain would trade safety for signal; leave on `stable` and fix warnings as they arise |
| `opal-sys` unexpectedly picks up `STAR/source/opal` on a self-hosted runner | All listed runners are GitHub-hosted; the tree won't exist. `build.rs` would compile it with `-mavx2` and likely fail on non-AVX2 runners — a nonissue today, documented here for future self-hosted runners |
| Draft release leaks if `publish-release` fails mid-run | Job is a single `gh release edit` call; rerun fixes it |

## File manifest (what the implementation will create/modify)

Created:
- `.github/workflows/ci.yml`
- `.github/workflows/release.yml`

Modified:
- `Cargo.toml` — update `repository` field only
