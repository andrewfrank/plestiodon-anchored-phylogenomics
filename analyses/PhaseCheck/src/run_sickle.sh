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

PROJECT_SUBDIR="sickle"
#$ -N sickle
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd
#$ -pe smp 8

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

DATA_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/concat_reads"

OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/trimmed_reads"
if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

# Run sickle
R1_TRIM="$DATAFILENAME"_R1_trim.fastq
R2_TRIM="$DATAFILENAME"_R2_trim.fastq
S_TRIM="$DATAFILENAME"_S_trim.fastq

sickle pe \
  -f $DATA_DIR/"$DATAFILENAME"_R1_cat.fastq.gz \
  -r $DATA_DIR/"$DATAFILENAME"_R2_cat.fastq.gz \
  -t sanger \
  -o $OUTPUT_DIR/$R1_TRIM \
  -p $OUTPUT_DIR/$R2_TRIM \
  -s $OUTPUT_DIR/$S_TRIM \
  -q 30 \
  -l 45
