#!/usr/bin/env python

# Planned script
# input
#   read ref fasta
#   read allele 1 fasta
#   read allele 2 fasta
# check if length(ref) = length(allele1) = length(allele2)
# run sequence comparisons
#   if ref[i] = allele1[i]      then c1[i] = 0, else c1[i] = 1
#   if ref[i] = allele2[i]      then c2[i] = 0, else c2[i] = 1
#   if allele1[i] = allele2[i]  then c3[i] = 0, else c3[i] = 1
# determine phase notation
#   if      c1[i] = 0 && c2[i] = 0,                 then phase = "0|0" & alt = ref[i]
#   elif    c1[i] = 0 && c2[i] = 1,                 then phase = "0|1" & alt = allele1[i]
#   elif    c1[i] = 1 && c2[i] = 0,                 then phase = "1|0" & alt = allele2[i]
#   elif    c1[i] = 1 && c2[i] = 1 && c3[i] = 0,    then phase = "1|1" & alt = allele1[i]
#   elif    c1[i] = 1 && c2[i] = 1 && c3[i] = 1,    then phase = "1|2" & alt = allele1[i],allele2[i]
# write to columns

from Bio import SeqIO
import os

input_path = "/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/data/original_seqs_byIndividual"
filename = "I7487_seq1.fasta"

full_path = os.path.join(input_path,filename)

# Maybe just read the whole set of 3 sequences into memory?

#for index, record in enumerate(SeqIO.parse(full_path,"fasta")):
#    print(record.seq)
