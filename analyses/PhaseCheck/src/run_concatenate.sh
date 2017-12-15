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

#$ -N concat
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

WORKING_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/original_reads"
cd "$WORKING_DIR"

# Concatenate read files
R1_CAT="$DATAFILENAME"_R1_cat.fastq.gz
R2_CAT="$DATAFILENAME"_R2_cat.fastq.gz

cat $WORKING_DIR/"$DATAFILENAME"_R1*.fastq.gz > "../concat_reads/$R1_CAT"
cat $WORKING_DIR/"$DATAFILENAME"_R2*.fastq.gz > "../concat_reads/$R2_CAT"
