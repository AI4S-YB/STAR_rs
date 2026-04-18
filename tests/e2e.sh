#!/usr/bin/env bash
# End-to-end regression tests for the Rust STAR port.
#
# For each covered code path, runs both the Rust binary and the C++ STAR
# reference on identical inputs and diffs selected outputs (SAM body,
# SJ.out.tab, sjdbInfo.txt, sjdbList.out.tab, genome index files).
#
# Usage:
#   ./tests/e2e.sh                        # run from star-rs/
#   STAR_REF=/path/STAR ./tests/e2e.sh    # override C++ STAR path
#   STAR_RS=target/release/star ./tests/e2e.sh
#   KEEP=1 ./tests/e2e.sh                 # keep scratch dir

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FIX="$SCRIPT_DIR/fixtures"

STAR_REF="${STAR_REF:-$ROOT/../STAR/bin/Linux_x86_64/STAR}"
STAR_RS="${STAR_RS:-$ROOT/target/release/star}"

SCRATCH="$(mktemp -d)"
trap '[ -n "${KEEP:-}" ] || rm -rf "$SCRATCH"' EXIT
[ -n "${KEEP:-}" ] && echo "Scratch: $SCRATCH"

FAIL=0
pass() { printf "  \033[32mPASS\033[0m %s\n" "$1"; }
fail() { printf "  \033[31mFAIL\033[0m %s\n" "$1"; FAIL=$((FAIL+1)); }

cmp_or() {
  local label="$1" a="$2" b="$3"
  if cmp -s "$a" "$b"; then pass "$label"; else fail "$label ($a vs $b)"; fi
}

cmp_sam_body() {
  local label="$1" a="$2" b="$3"
  if diff -q <(grep -v '^@' "$a") <(grep -v '^@' "$b") >/dev/null; then
    pass "$label"
  else
    fail "$label"
    diff <(grep -v '^@' "$a") <(grep -v '^@' "$b") | head -20
  fi
}

[ -x "$STAR_REF" ] || { echo "C++ STAR not found at $STAR_REF"; exit 2; }
[ -x "$STAR_RS" ] || { echo "Rust STAR not found at $STAR_RS (run: cargo build --release)"; exit 2; }

echo "== star-rs regression =="
echo "   rust: $STAR_RS"
echo "   c++ : $STAR_REF"

# ---------------------------------------------------------------------------
# Shared genome: tiny (4 KB chr1)
# ---------------------------------------------------------------------------
GEN_REF="$SCRATCH/gen_ref"
GEN_RS="$SCRATCH/gen_rs"
mkdir -p "$GEN_REF" "$GEN_RS"

echo
echo "[1] genomeGenerate bit-exact"
"$STAR_REF" --runMode genomeGenerate --runThreadN 1 \
  --genomeDir "$GEN_REF" --genomeFastaFiles "$FIX/tiny/chr1.fa" \
  --genomeSAindexNbases 6 --outFileNamePrefix "$GEN_REF/" >/dev/null 2>&1
"$STAR_RS" --runMode genomeGenerate --runThreadN 1 \
  --genomeDir "$GEN_RS" --genomeFastaFiles "$FIX/tiny/chr1.fa" \
  --genomeSAindexNbases 6 --outFileNamePrefix "$GEN_RS/" >/dev/null 2>&1

for f in Genome SA SAindex chrLength.txt chrName.txt chrNameLength.txt chrStart.txt; do
  cmp_or "genome/$f" "$GEN_RS/$f" "$GEN_REF/$f"
done

# ---------------------------------------------------------------------------
# alignReads SE, single-thread
# ---------------------------------------------------------------------------
echo
echo "[2] alignReads SE (runThreadN=1)"
R="$SCRATCH/se1_rs/"; C="$SCRATCH/se1_ref/"; mkdir -p "$R" "$C"
"$STAR_REF" --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outFileNamePrefix "$C" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outFileNamePrefix "$R" >/dev/null 2>&1
cmp_sam_body "se1/Aligned.out.sam body" "$R/Aligned.out.sam" "$C/Aligned.out.sam"
cmp_or       "se1/SJ.out.tab"           "$R/SJ.out.tab"      "$C/SJ.out.tab"

# ---------------------------------------------------------------------------
# alignReads PE, runThreadN=4
# ---------------------------------------------------------------------------
echo
echo "[3] alignReads PE (runThreadN=4)"
R="$SCRATCH/pe4_rs/"; C="$SCRATCH/pe4_ref/"; mkdir -p "$R" "$C"
"$STAR_REF" --runMode alignReads --runThreadN 4 \
  --genomeDir "$GEN_REF" \
  --readFilesIn "$FIX/tiny/reads_r1.fq" "$FIX/tiny/reads_r2.fq" \
  --outFileNamePrefix "$C" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 4 \
  --genomeDir "$GEN_REF" \
  --readFilesIn "$FIX/tiny/reads_r1.fq" "$FIX/tiny/reads_r2.fq" \
  --outFileNamePrefix "$R" >/dev/null 2>&1
cmp_sam_body "pe4/Aligned.out.sam body" "$R/Aligned.out.sam" "$C/Aligned.out.sam"
cmp_or       "pe4/SJ.out.tab"           "$R/SJ.out.tab"      "$C/SJ.out.tab"

# ---------------------------------------------------------------------------
# sjdbFileChrStartEnd (exercises sjdb_prepare + sjdb_build_index)
# ---------------------------------------------------------------------------
echo
echo "[4] alignReads --sjdbFileChrStartEnd (sjdb_build_index path)"
SJG_REF="$SCRATCH/sjg_ref"; SJG_RS="$SCRATCH/sjg_rs"; mkdir -p "$SJG_REF" "$SJG_RS"
"$STAR_REF" --runMode genomeGenerate --runThreadN 1 \
  --genomeDir "$SJG_REF" --genomeFastaFiles "$FIX/sjtest/chr.fa" \
  --genomeSAindexNbases 4 --outFileNamePrefix "$SJG_REF/" >/dev/null 2>&1
"$STAR_RS"  --runMode genomeGenerate --runThreadN 1 \
  --genomeDir "$SJG_RS"  --genomeFastaFiles "$FIX/sjtest/chr.fa" \
  --genomeSAindexNbases 4 --outFileNamePrefix "$SJG_RS/"  >/dev/null 2>&1
for f in Genome SA SAindex chrStart.txt; do
  cmp_or "sjtest-gen/$f" "$SJG_RS/$f" "$SJG_REF/$f"
done

R="$SCRATCH/sj_rs/"; C="$SCRATCH/sj_ref/"; mkdir -p "$R" "$C"
"$STAR_REF" --runMode alignReads --runThreadN 1 \
  --genomeDir "$SJG_REF" --readFilesIn "$FIX/sjtest/r.fq" \
  --sjdbFileChrStartEnd "$FIX/sjtest/sj.tab" --sjdbOverhang 30 \
  --outFileNamePrefix "$C" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 1 \
  --genomeDir "$SJG_REF" --readFilesIn "$FIX/sjtest/r.fq" \
  --sjdbFileChrStartEnd "$FIX/sjtest/sj.tab" --sjdbOverhang 30 \
  --outFileNamePrefix "$R" >/dev/null 2>&1
cmp_or "sjtest/_STARgenome/sjdbInfo.txt"     "$R/_STARgenome/sjdbInfo.txt"     "$C/_STARgenome/sjdbInfo.txt"
cmp_or "sjtest/_STARgenome/sjdbList.out.tab" "$R/_STARgenome/sjdbList.out.tab" "$C/_STARgenome/sjdbList.out.tab"

# ---------------------------------------------------------------------------
# --twopassMode Basic
# ---------------------------------------------------------------------------
echo
echo "[5] alignReads --twopassMode Basic (2-pass wiring)"
R="$SCRATCH/2p_rs/"; C="$SCRATCH/2p_ref/"; mkdir -p "$R" "$C"
"$STAR_REF" --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --twopassMode Basic --sjdbOverhang 50 \
  --outFileNamePrefix "$C" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --twopassMode Basic --sjdbOverhang 50 \
  --outFileNamePrefix "$R" >/dev/null 2>&1
cmp_sam_body "2p/Aligned.out.sam body"            "$R/Aligned.out.sam" "$C/Aligned.out.sam"
cmp_or       "2p/SJ.out.tab"                       "$R/SJ.out.tab"      "$C/SJ.out.tab"
cmp_or       "2p/_STARpass1/SJ.out.tab"            "$R/_STARpass1/SJ.out.tab" "$C/_STARpass1/SJ.out.tab"
cmp_or       "2p/_STARgenome/sjdbInfo.txt"         "$R/_STARgenome/sjdbInfo.txt"     "$C/_STARgenome/sjdbInfo.txt"
cmp_or       "2p/_STARgenome/sjdbList.out.tab"     "$R/_STARgenome/sjdbList.out.tab" "$C/_STARgenome/sjdbList.out.tab"

# ---------------------------------------------------------------------------
# Chimeric detection (--chimSegmentMin + --chimOutType Junctions)
# ---------------------------------------------------------------------------
echo
echo "[6] alignReads chimeric detection (legacy Old path)"
CHG_REF="$SCRATCH/chg_ref"; CHG_RS="$SCRATCH/chg_rs"; mkdir -p "$CHG_REF" "$CHG_RS"
"$STAR_REF" --runMode genomeGenerate --runThreadN 1 \
  --genomeDir "$CHG_REF" --genomeFastaFiles "$FIX/chimtest/chr.fa" \
  --genomeSAindexNbases 5 --outFileNamePrefix "$CHG_REF/" >/dev/null 2>&1

R="$SCRATCH/chim_rs/"; C="$SCRATCH/chim_ref/"; mkdir -p "$R" "$C"
"$STAR_REF" --runMode alignReads --runThreadN 1 \
  --genomeDir "$CHG_REF" --readFilesIn "$FIX/chimtest/r.fq" \
  --chimSegmentMin 20 --chimOutType Junctions \
  --outFileNamePrefix "$C" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 1 \
  --genomeDir "$CHG_REF" --readFilesIn "$FIX/chimtest/r.fq" \
  --chimSegmentMin 20 --chimOutType Junctions \
  --outFileNamePrefix "$R" >/dev/null 2>&1
cmp_sam_body "chim/Aligned.out.sam body"  "$R/Aligned.out.sam" "$C/Aligned.out.sam"
cmp_or       "chim/Chimeric.out.junction" "$R/Chimeric.out.junction" "$C/Chimeric.out.junction"

# runThreadN=4: deterministic chimeric output across worker chunks.
R4="$SCRATCH/chim_rs4/"; C4="$SCRATCH/chim_ref4/"; mkdir -p "$R4" "$C4"
"$STAR_REF" --runMode alignReads --runThreadN 4 \
  --genomeDir "$CHG_REF" --readFilesIn "$FIX/chimtest/r.fq" \
  --chimSegmentMin 20 --chimOutType Junctions \
  --outFileNamePrefix "$C4" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 4 \
  --genomeDir "$CHG_REF" --readFilesIn "$FIX/chimtest/r.fq" \
  --chimSegmentMin 20 --chimOutType Junctions \
  --outFileNamePrefix "$R4" >/dev/null 2>&1
cmp_or       "chim-mt/Chimeric.out.junction" "$R4/Chimeric.out.junction" "$C4/Chimeric.out.junction"

# ---------------------------------------------------------------------------
# --quantMode GeneCounts (M7.3)
# ---------------------------------------------------------------------------
echo
echo "[7] alignReads --quantMode GeneCounts"
GC_GEN="$SCRATCH/gc_gen"; mkdir -p "$GC_GEN"
# Genome with GTF must be built by C++ STAR (M7 does not yet port GTF
# parsing into --runMode genomeGenerate in Rust).
"$STAR_REF" --runMode genomeGenerate --runThreadN 1 \
  --genomeDir "$GC_GEN" --genomeFastaFiles "$FIX/gctest/chr.fa" \
  --genomeSAindexNbases 5 \
  --sjdbGTFfile "$FIX/gctest/ann.gtf" --sjdbOverhang 49 \
  --outFileNamePrefix "$GC_GEN/" >/dev/null 2>&1

R="$SCRATCH/gc_rs/"; C="$SCRATCH/gc_ref/"; mkdir -p "$R" "$C"
"$STAR_REF" --runMode alignReads --runThreadN 1 \
  --genomeDir "$GC_GEN" --readFilesIn "$FIX/gctest/r.fq" \
  --quantMode GeneCounts \
  --outFileNamePrefix "$C" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 1 \
  --genomeDir "$GC_GEN" --readFilesIn "$FIX/gctest/r.fq" \
  --quantMode GeneCounts \
  --outFileNamePrefix "$R" >/dev/null 2>&1
cmp_or "gc/ReadsPerGene.out.tab" "$R/ReadsPerGene.out.tab" "$C/ReadsPerGene.out.tab"

# runThreadN=4: deterministic GeneCounts across worker chunks.
R4="$SCRATCH/gc_rs4/"; C4="$SCRATCH/gc_ref4/"; mkdir -p "$R4" "$C4"
"$STAR_REF" --runMode alignReads --runThreadN 4 \
  --genomeDir "$GC_GEN" --readFilesIn "$FIX/gctest/r.fq" \
  --quantMode GeneCounts \
  --outFileNamePrefix "$C4" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 4 \
  --genomeDir "$GC_GEN" --readFilesIn "$FIX/gctest/r.fq" \
  --quantMode GeneCounts \
  --outFileNamePrefix "$R4" >/dev/null 2>&1
cmp_or "gc-mt/ReadsPerGene.out.tab" "$R4/ReadsPerGene.out.tab" "$C4/ReadsPerGene.out.tab"

# ---------------------------------------------------------------------------
# --outSAMtype BAM Unsorted (M7.4 minimum viable)
#
# Byte-exact BAM is out of scope for the minimum-viable port; we validate
# semantic equivalence instead:
#   * BAM is readable by samtools
#   * `samtools view` (body, sorted) matches C++ STAR for single and
#     multi-threaded runs
#   * no stray Aligned.out.sam.tmp remains
# ---------------------------------------------------------------------------
echo
echo "[8] alignReads --outSAMtype BAM Unsorted"
cmp_bam_body() {
  local label="$1" a="$2" b="$3"
  if ! command -v samtools >/dev/null 2>&1; then
    echo "  SKIP $label (samtools not available)"
    return 0
  fi
  local ta tb
  ta="$(mktemp)"; tb="$(mktemp)"
  samtools view "$a" | sort > "$ta"
  samtools view "$b" | sort > "$tb"
  if cmp -s "$ta" "$tb"; then pass "$label"; else fail "$label"; diff "$ta" "$tb" | head -10; fi
  rm -f "$ta" "$tb"
}

cmp_sorted_bam_body() {
  # Phase-2a acceptance: mapped records compared with within-coord qname
  # sort (algorithm parity without iRead coupling); unmapped records
  # compared as a multiset (sort | diff).
  local label="$1" a="$2" b="$3"
  if ! command -v samtools >/dev/null 2>&1; then
    echo "  SKIP $label (samtools not available)"
    return 0
  fi
  local am bm au bu
  am="$(mktemp)"; bm="$(mktemp)"; au="$(mktemp)"; bu="$(mktemp)"
  # Mapped: stable within-coord-bucket order by qname.
  samtools view -F 4 "$a" | awk 'BEGIN{FS=OFS="\t"} {print $3"_"$4, $1, $0}' | sort -k1,1 -k2,2 | cut -f3- > "$am"
  samtools view -F 4 "$b" | awk 'BEGIN{FS=OFS="\t"} {print $3"_"$4, $1, $0}' | sort -k1,1 -k2,2 | cut -f3- > "$bm"
  # Unmapped: multiset.
  samtools view -f 4 "$a" | sort > "$au"
  samtools view -f 4 "$b" | sort > "$bu"
  local ok=1
  if ! cmp -s "$am" "$bm"; then ok=0; echo "  DIFF mapped:"; diff "$am" "$bm" | head -8; fi
  if ! cmp -s "$au" "$bu"; then ok=0; echo "  DIFF unmapped:"; diff "$au" "$bu" | head -8; fi
  if [ "$ok" -eq 1 ]; then pass "$label"; else fail "$label"; fi
  rm -f "$am" "$bm" "$au" "$bu"
}

R="$SCRATCH/bam_rs/"; C="$SCRATCH/bam_ref/"; mkdir -p "$R" "$C"
"$STAR_REF" --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$C" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$R" >/dev/null 2>&1
cmp_bam_body "bam/Aligned.out.bam (body)" "$R/Aligned.out.bam" "$C/Aligned.out.bam"
# Scratch SAM must be deleted after conversion.
if [ ! -e "$R/Aligned.out.sam.tmp" ]; then pass "bam/no-scratch-sam"; else fail "bam/scratch SAM leaked"; fi

R4="$SCRATCH/bam_rs4/"; C4="$SCRATCH/bam_ref4/"; mkdir -p "$R4" "$C4"
"$STAR_REF" --runMode alignReads --runThreadN 4 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$C4" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 4 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$R4" >/dev/null 2>&1
cmp_bam_body "bam-mt/Aligned.out.bam (body)" "$R4/Aligned.out.bam" "$C4/Aligned.out.bam"

# ---------------------------------------------------------------------------
# --outSAMtype BAM SortedByCoordinate (phase 2a semantic match)
# ---------------------------------------------------------------------------
echo
echo "[9] alignReads --outSAMtype BAM SortedByCoordinate (phase 2a)"
R9="$SCRATCH/bamsort_rs/"; C9="$SCRATCH/bamsort_ref/"; mkdir -p "$R9" "$C9"
"$STAR_REF" --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$C9" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 1 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$R9" >/dev/null 2>&1
cmp_sorted_bam_body "bamsort/Aligned.sortedByCoord.out.bam (semantic)" \
  "$R9/Aligned.sortedByCoord.out.bam" "$C9/Aligned.sortedByCoord.out.bam"

R9m="$SCRATCH/bamsort_rs4/"; C9m="$SCRATCH/bamsort_ref4/"; mkdir -p "$R9m" "$C9m"
"$STAR_REF" --runMode alignReads --runThreadN 4 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$C9m" >/dev/null 2>&1
"$STAR_RS"  --runMode alignReads --runThreadN 4 \
  --genomeDir "$GEN_REF" --readFilesIn "$FIX/tiny/reads.fq" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$R9m" >/dev/null 2>&1
cmp_sorted_bam_body "bamsort-mt/Aligned.sortedByCoord.out.bam (semantic)" \
  "$R9m/Aligned.sortedByCoord.out.bam" "$C9m/Aligned.sortedByCoord.out.bam"

echo
if [ "$FAIL" -eq 0 ]; then
  echo "All regression tests passed."
  exit 0
else
  echo "$FAIL test(s) FAILED."
  exit 1
fi
