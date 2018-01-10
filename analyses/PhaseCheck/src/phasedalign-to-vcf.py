#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq
import os
import sys

sample =    sys.argv[1]
ref_file =  sys.argv[2]
al1_file =  sys.argv[3]
al2_file =  sys.argv[4]
out_file =  sys.argv[5]

base_dict = Seq.IUPAC.IUPACData.ambiguous_dna_values
std_bases = ["A","C","G","T"]
amb_bases = ["N","M","X","W","R","B","D","V","K","Y","S","H"]

refs = list(SeqIO.parse(ref_file, "fasta"))
al1s = list(SeqIO.parse(al1_file, "fasta"))
al2s = list(SeqIO.parse(al2_file, "fasta"))

header = ("##fileformat=VCFv4.2" + "\n"
            + "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" + "\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" + "\n"
            + "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another\">" + "\n")
colnames = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"
lines = [header,colnames]

for ref, al1, al2 in zip(refs, al1s, al2s):

    for i, (ref_base, al1_base, al2_base) in enumerate(zip(ref.seq, al1.seq, al2.seq)):

        chrom = ref.id
        pos = i + 1
        vcf_id = "."
        vcf_ref = ref_base
        qual = "."
        vcf_filter = "."
        info = "."

        if al1_base and al2_base not in amb_bases:
            vcf_format = "GT:PGT"
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

        if (al1_base in amb_bases) or (al2_base in amb_bases):
            vcf_format = "GT"
            dam_bases = base_dict.get(al1_base) + base_dict.get(al2_base)
            alt_bases = list(set(dam_bases))
            if len(alt_bases) == 2 and (ref_base in alt_bases):
                alt_bases.remove(ref_base)
                alt = ",".join(alt_bases)
                ind = "0/1"
            if len(alt_bases) == 2 and (ref_base not in alt_bases):
                alt = ",".join(alt_bases)
                ind = "1/2"
            if len(alt_bases) == 3 and (ref_base in alt_bases):
                alt_bases.remove(ref_base)
                alt = ",".join(alt_bases)
                ind = "0/1/2"
            if len(alt_bases) == 3 and (ref_base not in alt_bases):
                alt = ",".join(alt_bases)
                ind = "1/2/3"
            if len(alt_bases) == 4:
                alt_bases.remove(ref_base)
                alt = ",".join(alt_bases)
                ind = "0/1/2/3"

        line = chrom + "\t" + str(pos) + "\t" + vcf_id + "\t" + vcf_ref + "\t" + alt + "\t" + qual + "\t" + vcf_filter + "\t" + info + "\t" + vcf_format + "\t" + ind + "\n"
        lines.append(line)

vcf_out = open(out_file,"w")
vcf_out.writelines(lines)
vcf_out.close()
