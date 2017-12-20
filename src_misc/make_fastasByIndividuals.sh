#!/bin/bash

SAMPLES=(
I7487
I7488
I7489
I7490
I7491
I7492
I7493
I7495
I7496
I7497
I7498
I18066
I18068
I18069
)

ALLELES=(
seq1
seq2
)

INPUT_DIR="/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/data/original_seqs"
OUTPUT_DIR="/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/data/original_seqs_byIndividual"

cd $INPUT_DIR
awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' *.fasta > $OUTPUT_DIR/concat.fasta

cd $OUTPUT_DIR
for SAMPLE in "${SAMPLES[@]}"; do
  for ALLELE in "${ALLELES[@]}"; do
    grep -A 1 $SAMPLE concat.fasta |
    grep -A 1 $ALLELE |
    sed '/^--$/d' > "$SAMPLE"_"$ALLELE".fasta
  done
done
