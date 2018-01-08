#!/usr/bin/env python

from Bio import SeqIO
import os
import sys

sample =    sys.argv[1]
ref_file =  sys.argv[2]
al1_file =  sys.argv[3]
al2_file =  sys.argv[4]
out_file =  sys.argv[5]

refs = list(SeqIO.parse(ref_file, "fasta"))
al1s = list(SeqIO.parse(al1_file, "fasta"))
al2s = list(SeqIO.parse(al2_file, "fasta"))

header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"
lines = [header]

for ref, al1, al2 in zip(refs, al1s, al2s):

    for i, (ref_base, al1_base, al2_base) in enumerate(zip(ref.seq, al1.seq, al2.seq)):

        chrom = ref.id
        pos = i + 1
        vcf_id = "."
        vcf_ref = ref_base
        qual = "."
        vcf_filter = "."
        info = "."
        vcf_format = "GT,PGT"

        if ref_base == al1_base and ref_base == al2_base:
            alt = ref_base
            ind = "0/0:0|0"
        if ref_base == al1_base and ref_base != al2_base:
            alt = al2_base
            ind = "0/1:0|1"
        if ref_base != al1_base and ref_base == al2_base:
            alt = al1_base
            ind = "0/1:1|0"
        if ref_base != al1_base and ref_base != al2_base and al1_base == al2_base:
            alt = al1_base
            ind = "1/1:1|1"
        if ref_base != al1_base and ref_base != al2_base and al1_base != al2_base:
            alt = al1_base + "," + al2_base
            ind = "1/2:1|2"

        line = chrom + "\t" + str(pos) + "\t" + vcf_id + "\t" + vcf_ref + "\t" + alt + "\t" + qual + "\t" + vcf_filter + "\t" + info + "\t" + vcf_format + "\t" + ind + "\n"
        lines.append(line)

vcf_out = open(out_file,"w")
vcf_out.writelines(lines)
vcf_out.close()
