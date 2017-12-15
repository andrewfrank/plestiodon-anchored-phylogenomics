#!/bin/bash

DATA_FILES=(
I7487
#I7488
#I7489
#I7490
#I7491
#I7492
#I7493
#I7495
#I7496
#I7497
#I7498
#I18066
#I18068
#I18069
)

#$ -N callableloci
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

INDEX_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data/bwa-index"
DATA_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/filtered_mappings"
OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/coverage_filters"

if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

# Run GATK
module load GATK
java -jar $GATK37 \
  -T CallableLoci \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -I "$DATA_DIR"/"$DATAFILENAME"_fltrdmapping.bam \
  --minDepth 2 \
  -summary "$DATAFILENAME"_summary.txt \
  -o "$OUTPUT_DIR"/"$DATAFILENAME".bed
