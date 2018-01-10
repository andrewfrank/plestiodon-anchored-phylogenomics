#!/bin/bash

#SBATCH --job-name=UGT
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=andrew.frank@uconn.edu
#SBATCH -o UGT_%j.out
#SBATCH -e UGT_%j.err

module load GATK

cd "/home/CAM/afrank/plestiodon-anchored-phylogenomics/analyses/PhaseCheck"
java -jar $GATK \
 -T UnifiedGenotyper \
 -R ./data/bwa-index/ConsensusSeqs.fasta \
 -I ./data/bwa_mappings/I7487_mapping.bam \
 -I ./data/bwa_mappings/I7488_mapping.bam \
 -I ./data/bwa_mappings/I7489_mapping.bam \
 -I ./data/bwa_mappings/I7490_mapping.bam \
 -I ./data/bwa_mappings/I7491_mapping.bam \
 -I ./data/bwa_mappings/I7492_mapping.bam \
 -I ./data/bwa_mappings/I7493_mapping.bam \
 -I ./data/bwa_mappings/I7495_mapping.bam \
 -I ./data/bwa_mappings/I7496_mapping.bam \
 -I ./data/bwa_mappings/I7497_mapping.bam \
 -I ./data/bwa_mappings/I7498_mapping.bam \
 -I ./data/bwa_mappings/I18066_mapping.bam \
 -I ./data/bwa_mappings/I18068_mapping.bam \
 -I ./data/bwa_mappings/I18069_mapping.bam \
 -o ./data/variants_UGC/VariantCalls_frank-UGT.vcf
