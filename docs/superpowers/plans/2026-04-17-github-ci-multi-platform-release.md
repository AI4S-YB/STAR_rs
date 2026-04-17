# GitHub CI & Multi-Platform Release Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add GitHub Actions workflows that run fmt/clippy/test/MSRV on every PR and push, and build multi-platform release binaries on version tags or manual dispatch.

**Architecture:** Two workflow files under `.github/workflows/`. `ci.yml` runs four parallel jobs on push/PR to `main`. `release.yml` creates a draft GitHub Release, builds binaries across a 5-target matrix, and publishes the release only after all targets succeed.

**Tech Stack:** GitHub Actions, `dtolnay/rust-toolchain`, `Swatinem/rust-cache@v2`, `taiki-e/upload-rust-binary-action@v1`, `gh` CLI.

**Design spec:** `docs/superpowers/specs/2026-04-17-github-ci-multi-platform-release-design.md`

---

## File Structure

Files created:
- `.github/workflows/ci.yml` — CI checks (fmt, clippy, test matrix, MSRV)
- `.github/workflows/release.yml` — multi-platform release build + publish

File modified:
- `Cargo.toml` — update `repository` field

---

## Validation approach

GitHub Actions workflows are validated in three layers:

1. **YAML syntax** — run `python3 -c "import yaml; yaml.safe_load(open('<file>'))"` to catch parse errors. This runs locally without extra tooling.
2. **Action schema** — optional: install `actionlint` via `bash <(curl -fsSL https://raw.githubusercontent.com/rhysd/actionlint/main/scripts/download-actionlint.bash)` for richer linting. Not required; skip if `curl` is unavailable.
3. **Runtime** — the only way to fully verify a workflow is to run it on GitHub. That happens after merge: CI is exercised by the PR itself, and Release is exercised by pushing a preview tag.

Each task uses layer 1 as a hard gate. Layer 2 is a soft recommendation. Layer 3 is acknowledged but out of scope for the implementation session.

---

## Task 1: Update `Cargo.toml` repository field

**Files:**
- Modify: `Cargo.toml:23`

- [ ] **Step 1: Make the edit**

Change line 23 from:
```toml
repository = "https://github.com/example/star-rs"
```
to:
```toml
repository = "https://github.com/AI4S-YB/STAR_rs"
```

- [ ] **Step 2: Verify the workspace still parses**

Run: `cargo metadata --no-deps --format-version 1 >/dev/null`
Expected: exits 0 with no output.

- [ ] **Step 3: Commit**

```bash
git add Cargo.toml
git commit -m "chore: point repository field at AI4S-YB/STAR_rs"
```

---

## Task 2: Create `.github/workflows/ci.yml`

**Files:**
- Create: `.github/workflows/ci.yml`

- [ ] **Step 1: Create the directory**

Run: `mkdir -p .github/workflows`

- [ ] **Step 2: Write the workflow**

Create `.github/workflows/ci.yml` with this exact content:

```yaml
name: CI

on:
  pull_request:
    branches: [main]
  push:
    branches: [main]

concurrency:
  group: ci-${{ github.ref }}
  cancel-in-progress: true

env:
  CARGO_TERM_COLOR: always

jobs:
  fmt:
    name: Format
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: rustfmt
      - run: cargo fmt --all -- --check

  clippy:
    name: Clippy
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: clippy
      - uses: Swatinem/rust-cache@v2
        with:
          key: clippy
      - run: cargo clippy --workspace --all-targets -- -D warnings

  test:
    name: Test (${{ matrix.os }})
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - uses: Swatinem/rust-cache@v2
        with:
          key: test-${{ matrix.os }}
      - run: cargo test --workspace

  msrv:
    name: MSRV (1.75)
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@1.75.0
      - uses: Swatinem/rust-cache@v2
        with:
          key: msrv
      - run: cargo check --workspace --all-targets
```

- [ ] **Step 3: Validate YAML syntax**

Run: `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))"`
Expected: exits 0, no output. If Python 3 is unavailable, try `python -c ...`.

- [ ] **Step 4 (optional): Run actionlint**

If `actionlint` is on PATH, run `actionlint .github/workflows/ci.yml` and expect zero output. Otherwise install it with:
```bash
bash <(curl -fsSL https://raw.githubusercontent.com/rhysd/actionlint/main/scripts/download-actionlint.bash)
./actionlint .github/workflows/ci.yml
rm actionlint
```
Expected: no output. If install fails (no network, no curl), skip this step — YAML validity from Step 3 is sufficient to proceed.

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: add fmt/clippy/test/MSRV workflow"
```

---

## Task 3: Create `.github/workflows/release.yml`

**Files:**
- Create: `.github/workflows/release.yml`

- [ ] **Step 1: Write the workflow**

Create `.github/workflows/release.yml` with this exact content:

```yaml
name: Release

on:
  push:
    tags: ['v*']
  workflow_dispatch:
    inputs:
      tag:
        description: 'Tag to build (e.g. v0.0.0-dev)'
        required: true
        default: 'v0.0.0-dev'

permissions:
  contents: write

env:
  CARGO_TERM_COLOR: always

jobs:
  create-release:
    name: Create draft release
    runs-on: ubuntu-latest
    outputs:
      tag: ${{ steps.tag.outputs.tag }}
    steps:
      - uses: actions/checkout@v4
      - name: Resolve tag
        id: tag
        run: |
          if [ "${{ github.event_name }}" = "workflow_dispatch" ]; then
            TAG="${{ inputs.tag }}"
          else
            TAG="${GITHUB_REF#refs/tags/}"
          fi
          echo "tag=$TAG" >> "$GITHUB_OUTPUT"
      - name: Create draft release if missing
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          TAG: ${{ steps.tag.outputs.tag }}
        run: |
          if gh release view "$TAG" >/dev/null 2>&1; then
            echo "Release $TAG already exists, skipping create"
          else
            gh release create "$TAG" --draft --generate-notes --title "$TAG"
          fi

  build:
    name: Build (${{ matrix.target }})
    needs: create-release
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
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
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          targets: ${{ matrix.target }}
      - uses: Swatinem/rust-cache@v2
        with:
          key: release-${{ matrix.target }}
      - name: Install apt deps
        if: matrix.apt != ''
        run: |
          sudo apt-get update
          sudo apt-get install -y ${{ matrix.apt }}
      - uses: taiki-e/upload-rust-binary-action@v1
        with:
          bin: star
          target: ${{ matrix.target }}
          archive: star-${{ needs.create-release.outputs.tag }}-${{ matrix.target }}
          checksum: sha256
          tar: unix
          zip: windows
          ref: refs/tags/${{ needs.create-release.outputs.tag }}
          token: ${{ secrets.GITHUB_TOKEN }}

  publish-release:
    name: Publish release
    needs: [create-release, build]
    runs-on: ubuntu-latest
    steps:
      - name: Publish
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          TAG: ${{ needs.create-release.outputs.tag }}
        run: gh release edit "$TAG" --draft=false
```

- [ ] **Step 2: Validate YAML syntax**

Run: `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/release.yml'))"`
Expected: exits 0, no output.

- [ ] **Step 3 (optional): Run actionlint**

If `actionlint` is on PATH (or installable per Task 2 Step 4), run:
```bash
actionlint .github/workflows/release.yml
```
Expected: no output. Skip silently if not installable.

- [ ] **Step 4: Sanity-check matrix coverage**

Confirm the workflow file contains all five target triples. Run:
```bash
grep -E 'target: (x86_64-unknown-linux-gnu|x86_64-unknown-linux-musl|x86_64-apple-darwin|aarch64-apple-darwin|x86_64-pc-windows-msvc)$' .github/workflows/release.yml | wc -l
```
Expected: `5`.

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/release.yml
git commit -m "ci: add multi-platform release workflow"
```

---

## Task 4: Final verification

**Files:** none modified.

- [ ] **Step 1: Inspect the tree**

Run: `ls -la .github/workflows/`
Expected output includes `ci.yml` and `release.yml`.

- [ ] **Step 2: Confirm no other changes**

Run: `git status`
Expected: `nothing to commit, working tree clean` (assuming all earlier commits succeeded).

- [ ] **Step 3: Review recent commits**

Run: `git log --oneline -5`
Expected to show (most recent first):
1. `ci: add multi-platform release workflow`
2. `ci: add fmt/clippy/test/MSRV workflow`
3. `chore: point repository field at AI4S-YB/STAR_rs`
4. `docs: add design spec for GitHub CI & multi-platform release`
5. (prior repo commit)

- [ ] **Step 4: Summary for user**

Print a short note to the user including:
- How to exercise the workflows (open a PR for CI; push `v0.1.0` or use workflow_dispatch for Release).
- Reminder that runtime behavior can only be confirmed on GitHub.
- The expected six release assets (5 archives + a Windows `.zip`, each with a `.sha256` sidecar).

---

## Out of scope (explicit)

- README badges — deferred until workflows are green on GitHub.
- Publishing to crates.io.
- Docker images, code coverage, benchmarks.
- Self-hosted runners (would re-enable the AVX2 build path in `opal-sys`; flagged in the spec).
