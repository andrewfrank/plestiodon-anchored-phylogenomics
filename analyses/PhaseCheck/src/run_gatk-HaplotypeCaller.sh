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

#$ -N haplotypecaller
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

INDEX_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data/bwa-index"
DATA_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/filtered_mappings"
OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/gatk-HaplotypeCaller_gvcfs"

if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

# Run GATK
module load GATK
java -jar $GATK37 \
  -T HaplotypeCaller \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -I "$DATA_DIR"/"$DATAFILENAME"_fltrdmapping.bam \
  --emitRefConfidence GVCF \
  -o "$OUTPUT_DIR"/"$DATAFILENAME".g.vcf \
  -bamout "$OUTPUT_DIR"/"$DATAFILENAME"_haplotypes.bam
