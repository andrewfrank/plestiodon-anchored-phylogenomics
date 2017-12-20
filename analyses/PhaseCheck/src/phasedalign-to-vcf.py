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
#   if      c1[i] = 0 && c2[i] = 0,                 then phase = "0|0"
#   elif    c1[i] = 0 && c2[i] = 1,                 then phase = "0|1"
#   elif    c1[i] = 1 && c2[i] = 0,                 then phase = "1|0"
#   elif    c1[i] = 1 && c2[i] = 1 && c3[i] = 0,    then phase = "1|1"
#   elif    c1[i] = 1 && c2[i] = 1 && c3[i] = 1,    then phase = "1|2"
# determine alt allele notation
