#!/usr/bin/env python3
"""Build a tiny 2-chromosome genome + chimeric read fixture.

The chimeric read is chr1[100..150) + chr2[50..100) — a 100-nt fusion that
straddles a cross-chromosome junction. Running STAR with
`--chimSegmentMin 20 --chimOutType Junctions` should emit one line in
Chimeric.out.junction describing the breakpoint between chr1:150(+) and
chr2:50(+).
"""
import os
import random

random.seed(17)
OUT = os.path.dirname(os.path.abspath(__file__))


def randseq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


chr1 = randseq(1000)
chr2 = randseq(1000)

with open(os.path.join(OUT, "chr.fa"), "w") as f:
    f.write(">chr1\n")
    for i in range(0, len(chr1), 60):
        f.write(chr1[i:i + 60] + "\n")
    f.write(">chr2\n")
    for i in range(0, len(chr2), 60):
        f.write(chr2[i:i + 60] + "\n")

reads = []
# Exonic controls (non-chimeric, fall entirely within chr1 or chr2).
reads.append(("ex_chr1_0", chr1[200:260]))
reads.append(("ex_chr2_0", chr2[200:260]))
# Chimeric: fuse chr1 and chr2 at mid-span. 50 nt from each side.
reads.append(("chim_0", chr1[500:550] + chr2[500:550]))
reads.append(("chim_1", chr1[600:660] + chr2[600:660]))

with open(os.path.join(OUT, "r.fq"), "w") as f:
    for name, seq in reads:
        f.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
