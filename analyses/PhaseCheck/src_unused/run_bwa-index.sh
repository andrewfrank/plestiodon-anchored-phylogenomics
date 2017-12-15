#!/bin/bash

DATA_FILES=(
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

#$ -N bwa-mem
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

# Specify directories
WORKING_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data/bwa-index"

# Create output directory if not present
if [ ! -d "$WORKING_DIR" ]; then
mkdir -p "$WORKING_DIR"
fi

# Run bowtie2
bwa index \
  -p "$WORKING_DIR"/"$DATAFILENAME"_seqs \
  -a is \
  "$WORKING_DIR"/"$DATAFILENAME"_seqs.fasta
