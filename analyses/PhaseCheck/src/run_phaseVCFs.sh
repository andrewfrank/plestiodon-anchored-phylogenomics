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

#$ -N whatshap
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

INDEX_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data/bwa-index"
BAM_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/filtered_mappings"
VCF_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/variants"
OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/phased_vcfs"

if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

# Run WhatsHap
whatshap phase \
  --reference "$INDEX_DIR"/ConsensusSeqs.fasta \
  "$VCF_DIR"/"$DATAFILENAME".vcf \
  "$BAM_DIR"/"$DATAFILENAME"_fltrdmapping.bam \
  --indels \
  --ignore-read-groups \
  -o "$OUTPUT_DIR"/"$DATAFILENAME"_whp.vcf

# Make phased fastas
bgzip phased.vcf
tabix phased.vcf.gz
bcftools consensus -H 1 -f reference.fasta phased.vcf.gz > haplotype1.fasta
bcftools consensus -H 2 -f reference.fasta phased.vcf.gz > haplotype2.fasta
