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

#$ -N selectvariants
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

INDEX_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data/bwa-index"
WORKING_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/variants"

if [ ! -d "$WORKING_DIR" ]; then
mkdir -p "$WORKING_DIR"
fi

# Run GATK
module load GATK
java -jar $GATK37 \
  -T SelectVariants \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$WORKING_DIR"/VariantCalls.vcf \
  -o "$WORKING_DIR"/"$DATAFILENAME".vcf \
  --excludeFiltered \
  --excludeNonVariants \
  --removeUnusedAlternates \
  -sn "$DATAFILENAME"

bgzip -c \
  "$WORKING_DIR"/"$DATAFILENAME".vcf \
  > "$WORKING_DIR"/"$DATAFILENAME".vcf.gz
tabix -p vcf "$WORKING_DIR"/"$DATAFILENAME".vcf.gz

java -jar $GATK37 \
  -T SelectVariants \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$WORKING_DIR"/VariantCalls_fltrd.vcf \
  -o "$WORKING_DIR"/"$DATAFILENAME"_fltrd.vcf \
  --excludeFiltered \
  --excludeNonVariants \
  --removeUnusedAlternates \
  -sn "$DATAFILENAME"

bgzip -c \
  "$WORKING_DIR"/"$DATAFILENAME"_fltrd.vcf \
  > "$WORKING_DIR"/"$DATAFILENAME"_fltrd.vcf.gz
tabix -p vcf "$WORKING_DIR"/"$DATAFILENAME"_fltrd.vcf.gz
