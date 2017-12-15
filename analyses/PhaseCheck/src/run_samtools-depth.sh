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

#$ -N samtools
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"
DATA_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/picard_grpdmappings"
OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/samtools_coverage"

if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

# Run samtools
samtools depth \
  $DATA_DIR/"$DATAFILENAME"_grpdmapping.bam \
  > $OUTPUT_DIR/"$DATAFILENAME"_grpdmapping.coverage
