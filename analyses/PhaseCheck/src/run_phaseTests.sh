#!/bin/bash

#$ -N phasetest
#$ -t 1-1
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME=I7487

INDEX_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data/bwa-index"
BAM_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/filtered_mappings"
VCF_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/variants"
OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/phasing_tests"

if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

if [ ! -d "$OUTPUT_DIR"/"$DATAFILENAME"_psr ]; then
mkdir -p "$OUTPUT_DIR"/"$DATAFILENAME"_psr
fi

# Run phASER
#   (--gw_phase_vcf 2) puts phase information into the GT field of the vcf
#   (--as_q_cutoff 0) solves a bug with the AS tag - I should go make sure my
#     BAM files' AS tags are OK
#   (--pass_only 0) allows my variants that weren't filtered to go through
python ~/develop/phaser/phaser/phaser.py \
  --bam "$BAM_DIR"/"$DATAFILENAME"_fltrdmapping.bam \
  --vcf "$VCF_DIR"/"$DATAFILENAME".vcf.gz \
  --sample "$DATAFILENAME" \
  --mapq 60 \
  --baseq 10 \
  --paired_end 1 \
  --as_q_cutoff 0 \
  --pass_only 0 \
  --include_indels 1 \
  --id_separator "-" \
  --python_string python2.7 \
  --gw_phase_vcf 2 \
  --o "$OUTPUT_DIR"/"$DATAFILENAME"_psr/"$DATAFILENAME"_psr

# Run HapCUT2
~/develop/HapCUT2/build/extractHAIRS \
  --bam "$BAM_DIR"/"$DATAFILENAME"_fltrdmapping.bam \
  --VCF "$VCF_DIR"/"$DATAFILENAME".vcf \
  --ref "$INDEX_DIR"/ConsensusSeqs.fasta \
  --indels 1 \
  --out "$OUTPUT_DIR"/"$DATAFILENAME"_hc2.frag
~/develop/HapCUT2/build/HAPCUT2 \
  --fragments "$OUTPUT_DIR"/"$DATAFILENAME"_hc2.frag \
  --VCF "$VCF_DIR"/"$DATAFILENAME".vcf \
  --out "$OUTPUT_DIR"/"$DATAFILENAME"_hc2.out
# At this point I converted hc2.out to a vcf using fgbio on my own computer

# Run HapCompass
~/develop/jdk1.8.0_91/jre/bin/java \
  -jar ~/develop/hapcompass_v0.8.2/hapcompass.jar \
  -Xmx1g \
  -bam "$BAM_DIR"/"$DATAFILENAME"_fltrdmapping.bam \
  -vcf "$VCF_DIR"/"$DATAFILENAME".vcf \
  -o "$OUTPUT_DIR"/"$DATAFILENAME"_hcp/"$DATAFILENAME"
~/develop/jdk1.8.0_91/jre/bin/java \
  -jar ~/develop/hapcompass_v0.8.2/hc2vcf.jar \
  "$OUTPUT_DIR"/"$DATAFILENAME"_hcp/"$DATAFILENAME"_MWER_solution.txt \
  "$VCF_DIR"/"$DATAFILENAME"_hcp.vcf \
  2 \
  true

# Run WhatsHap
whatshap phase \
  --reference "$INDEX_DIR"/ConsensusSeqs.fasta \
  "$VCF_DIR"/"$DATAFILENAME".vcf \
  "$BAM_DIR"/"$DATAFILENAME"_fltrdmapping.bam \
  --indels \
  --ignore-read-groups \
  -o "$OUTPUT_DIR"/"$DATAFILENAME"_whp.vcf

## phASER: PHASED  587 of 881 all variants (= 0.666288) with at least one other variant
## WhatsHap: Phased:      667 (     605 SNVs)
