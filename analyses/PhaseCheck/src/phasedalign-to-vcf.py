#!/usr/bin/env python

# Planned script
# input
#   read ref fasta
#   read allele 1 fasta
#   read allele 2 fasta
# check if length(ref) = length(allele1) = length(allele2)
# run sequence comparisons
#   if ref[j] = allele1[j]      then c1[j] = 0, else c1[j] = 1
#   if ref[j] = allele2[j]      then c2[j] = 0, else c2[j] = 1
#   if allele1[j] = allele2[j]  then c3[j] = 0, else c3[j] = 1
# determine phase notation
#   if      c1[j] = 0 && c2[j] = 0,                 then phase = "0|0" & alt = ref[j]
#   elif    c1[j] = 0 && c2[j] = 1,                 then phase = "0|1" & alt = allele2[j]
#   elif    c1[j] = 1 && c2[j] = 0,                 then phase = "1|0" & alt = allele1[j]
#   elif    c1[j] = 1 && c2[j] = 1 && c3[j] = 0,    then phase = "1|1" & alt = allele1[j]
#   elif    c1[j] = 1 && c2[j] = 1 && c3[j] = 1,    then phase = "1|2" & alt = allele1[j],allele2[j]
# write to columns

from Bio import SeqIO
import os

working_dir = "/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/data"

sample = "I7487"

ref_file = os.path.join(working_dir,"consensus_seqs","ConsensusSeqs_noAmbi.fasta")
al1_file = os.path.join(working_dir,"original_seqs_byIndividual",sample + "_seq1.fasta")
al2_file = os.path.join(working_dir,"original_seqs_byIndividual",sample + "_seq2.fasta")
out_file = os.path.join(working_dir,"original_vcfs",sample + "_original.vcf")

ref = list(SeqIO.parse(ref_file, "fasta"))
al1 = list(SeqIO.parse(al1_file, "fasta"))
al2 = list(SeqIO.parse(al2_file, "fasta"))

chrom = []
pos = []
vcf_id = []
vcf_ref = []
alt = []
qual = []
vcf_filter = []
info = []
vcf_format = []
ind = []

# Start for loop here over each locus

for i in range( len(ref) - 1):

    for j in range( len(ref[i].seq) - 1 ):
        chrom.append(ref[i].id)
        pos.append(j + 1)
        vcf_id.append(".")
        vcf_ref.append(ref[i].seq[j])
        qual.append(".")
        vcf_filter.append(".")
        info.append(".")
        vcf_format.append("GT,PGT")

    # ternary operators:
    # https://www.pythoncentral.io/one-line-if-statement-in-python-ternary-conditional-operator/

    c1 = []
    c2 = []
    c3 = []

    for j in range( len(ref[i].seq) - 1 ):
        c1.append(0 if ref[i].seq[j] == al1[i].seq[j] else 1)
        c2.append(0 if ref[i].seq[j] == al2[i].seq[j] else 1)
        c3.append(0 if al1[i].seq[j] == al2[i].seq[j] else 1)

    for j in range( len(ref[i].seq) - 1 ):
        if c1[j] == 0 and c2[j] == 0:
            alt.append(ref[i].seq[j])
            ind.append("0/0:0|0")
        if c1[j] == 0 and c2[j] == 1:
            alt.append(al2[i].seq[j])
            ind.append("0/1:0|1")
        if c1[j] == 1 and c2[j] == 0:
            alt.append(al1[i].seq[j])
            ind.append("0/1:1|0")
        if c1[j] == 1 and c2[j] == 1 and c3[j] == 0:
            alt.append(al1[i].seq[j])
            ind.append("1/1:1|1")
        if c1[j] == 1 and c2[j] == 1 and c3[j] == 1:
            alt.append(al1[i].seq[j] + "," + al2[i].seq[j])
            ind.append("1/2:1|2")

# Separate loop here
lines = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"]
for i in range( len(chrom) - 1):
    lines.append(chrom[i] + "\t" +
                str(pos[i]) + "\t" +
                vcf_id[i] + "\t" +
                vcf_ref[i] + "\t" +
                alt[i] + "\t" +
                qual[i] + "\t" +
                vcf_filter[i] + "\t" +
                info[i] + "\t" +
                vcf_format[i] + "\t" +
                ind[i] + "\n")

vcf_out = open(out_file,"w")
vcf_out.writelines(lines)
vcf_out.close()
