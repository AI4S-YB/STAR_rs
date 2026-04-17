#!/usr/bin/env python3
"""Build a tiny genome + GTF + reads to exercise --quantMode GeneCounts.

Layout:
  chr1 (2000 nt)
    geneA: exon [100..199] (+)
           exon [300..399] (+)
    geneB: exon [800..899] (-)
    geneC: exon [1200..1299] (+)
           exon [1400..1499] (+)

Reads (all 50 nt, all length 50):
  - ex_A:   chr1[110..160)   -> inside geneA exon1 (strand +)
  - ex_A2:  chr1[310..360)   -> inside geneA exon2 (strand +)
  - ex_B:   chr1[810..860)   -> inside geneB exon1 (strand -)
  - ex_C:   chr1[1210..1260) -> inside geneC exon1 (strand +)
  - ig_1:   chr1[500..550)   -> intergenic (no feature)
  - ig_2:   chr1[1000..1050) -> intergenic (no feature)
"""
import os
import random

random.seed(11)
OUT = os.path.dirname(os.path.abspath(__file__))


def randseq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


seq = randseq(2000)
with open(os.path.join(OUT, "chr.fa"), "w") as f:
    f.write(">chr1\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i + 60] + "\n")


def gtf_line(feature, start, end, strand, attrs):
    return "\t".join([
        "chr1", "mk", feature,
        str(start), str(end), ".", strand, ".",
        attrs,
    ]) + "\n"


with open(os.path.join(OUT, "ann.gtf"), "w") as f:
    # geneA: two exons on +, transcript tA
    gene_attrs = 'gene_id "geneA"; gene_name "geneA_n"; gene_biotype "protein_coding";'
    tr_attrs = 'gene_id "geneA"; transcript_id "tA"; gene_name "geneA_n"; gene_biotype "protein_coding";'
    f.write(gtf_line("gene", 100, 400, "+", gene_attrs))
    f.write(gtf_line("transcript", 100, 400, "+", tr_attrs))
    f.write(gtf_line("exon", 100, 199, "+", tr_attrs + ' exon_number "1";'))
    f.write(gtf_line("exon", 300, 399, "+", tr_attrs + ' exon_number "2";'))

    # geneB: one exon on -
    gene_attrs = 'gene_id "geneB"; gene_name "geneB_n"; gene_biotype "protein_coding";'
    tr_attrs = 'gene_id "geneB"; transcript_id "tB"; gene_name "geneB_n"; gene_biotype "protein_coding";'
    f.write(gtf_line("gene", 800, 900, "-", gene_attrs))
    f.write(gtf_line("transcript", 800, 900, "-", tr_attrs))
    f.write(gtf_line("exon", 800, 899, "-", tr_attrs + ' exon_number "1";'))

    # geneC: two exons on +
    gene_attrs = 'gene_id "geneC"; gene_name "geneC_n"; gene_biotype "protein_coding";'
    tr_attrs = 'gene_id "geneC"; transcript_id "tC"; gene_name "geneC_n"; gene_biotype "protein_coding";'
    f.write(gtf_line("gene", 1200, 1500, "+", gene_attrs))
    f.write(gtf_line("transcript", 1200, 1500, "+", tr_attrs))
    f.write(gtf_line("exon", 1200, 1299, "+", tr_attrs + ' exon_number "1";'))
    f.write(gtf_line("exon", 1400, 1499, "+", tr_attrs + ' exon_number "2";'))


reads = []
# 1-based GTF, 0-based python: exon 100..199 (incl) == seq[99:199]
reads.append(("ex_A_1",  seq[110:160]))
reads.append(("ex_A_2",  seq[310:360]))
reads.append(("ex_A_3",  seq[120:170]))
reads.append(("ex_B_1",  seq[810:860]))
reads.append(("ex_B_2",  seq[820:870]))
reads.append(("ex_C_1",  seq[1210:1260]))
reads.append(("ig_1",    seq[500:550]))
reads.append(("ig_2",    seq[1000:1050]))

with open(os.path.join(OUT, "r.fq"), "w") as f:
    for name, s in reads:
        f.write(f"@{name}\n{s}\n+\n{'I' * len(s)}\n")
