# Build a tiny 2-exon "gene" with canonical GT..AG intron and
# generate reads spanning the junction.
import random, os
random.seed(1)

# Base sequence from chr1: use random bases so the junction is easy to
# detect via AG/GT motif (force it).
def randseq(n):
    return "".join(random.choice("ACGT") for _ in range(n))

exon1 = randseq(200)
intron = "GT" + randseq(196) + "AG"   # 200 bp
exon2 = randseq(200)
ref = exon1 + intron + exon2 + randseq(600)   # pad to 1200

with open("/tmp/star_sjtest/chr.fa","w") as f:
    f.write(">chr1\n")
    for i in range(0, len(ref), 60):
        f.write(ref[i:i+60]+"\n")

# Reads: spanning junction - 50mer with 25nt from end of exon1 and 25nt from start of exon2
with open("/tmp/star_sjtest/r.fq","w") as f:
    for i,off in enumerate([175,180,185,190,195]):
        seq = exon1[off:] + exon2[:50 - (200-off)]
        f.write(f"@sp_{i}\n{seq}\n+\n{'I'*len(seq)}\n")
    # Also normal exonic reads
    for i,off in enumerate([0, 50, 400, 500]):
        seq = ref[off:off+50]
        f.write(f"@ex_{i}\n{seq}\n+\n{'I'*len(seq)}\n")
