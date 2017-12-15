#!/bin/bash

#$ -N filtervariants
#$ -t 1-1
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

INDEX_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data/bwa-index"
DATA_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/gatk-HaplotypeCaller_gvcfs"
OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/variants"

if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

# Run GATK
module load GATK
java -jar $GATK37 \
  -T GenotypeGVCFs \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  --variant "$DATA_DIR"/I7487.g.vcf \
  --variant "$DATA_DIR"/I7488.g.vcf \
  --variant "$DATA_DIR"/I7489.g.vcf \
  --variant "$DATA_DIR"/I7490.g.vcf \
  --variant "$DATA_DIR"/I7491.g.vcf \
  --variant "$DATA_DIR"/I7492.g.vcf \
  --variant "$DATA_DIR"/I7493.g.vcf \
  --variant "$DATA_DIR"/I7495.g.vcf \
  --variant "$DATA_DIR"/I7496.g.vcf \
  --variant "$DATA_DIR"/I7497.g.vcf \
  --variant "$DATA_DIR"/I7498.g.vcf \
  --variant "$DATA_DIR"/I18066.g.vcf \
  --variant "$DATA_DIR"/I18068.g.vcf \
  --variant "$DATA_DIR"/I18069.g.vcf \
  -o "$OUTPUT_DIR"/VariantCalls.vcf

java -jar $GATK37 \
  -T SelectVariants \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$OUTPUT_DIR"/VariantCalls.vcf \
  -selectType SNP \
  -o "$OUTPUT_DIR"/VariantCalls_SNPsOnly.vcf

# I changed the filter expression because it was throwing an error on some
# variants where one of the criteria didn't exist, and thus (I think) threw output
# the entire SNP despite it passing other criteria
# According to this thread: http://seqanswers.com/forums/showthread.php?t=16556
# you can input the filter expressions individually and it should evaluate
# properly
# Here's the original GATK suggested filter expression:
#   --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
#   --filterName "gatk_default_SNP_filter" \
java -jar $GATK37 \
  -T VariantFiltration \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$OUTPUT_DIR"/VariantCalls_SNPsOnly.vcf \
  --filterExpression "QD < 2.0" \
  --filterName "QDfilter" \
  --filterExpression "FS > 60.0" \
  --filterName "FSfilter" \
  --filterExpression "MQ < 40.0" \
  --filterName "MQfilter" \
  --filterExpression "MQRankSum < -12.5" \
  --filterName "MQRankSumfilter" \
  --filterExpression "ReadPosRankSum < -8.0" \
  --filterName "ReadPosRankSumfilter" \
  -o "$OUTPUT_DIR"/VariantCalls_fltrdSNPs.vcf

java -jar $GATK37 \
  -T SelectVariants \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$OUTPUT_DIR"/VariantCalls.vcf \
  -selectType INDEL \
  -o "$OUTPUT_DIR"/VariantCalls_INDELsOnly.vcf

# Here's the original GATK suggested filter expression:
#   --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
#   --filterName "gatk_default_SNP_filter" \
java -jar $GATK37 \
  -T VariantFiltration \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$OUTPUT_DIR"/VariantCalls_INDELsOnly.vcf \
  --filterExpression "QD < 2.0" \
  --filterName "QDfilter" \
  --filterExpression "FS > 200.0" \
  --filterName "FSfilter" \
  --filterExpression "ReadPosRankSum < -20.0" \
  --filterName "ReadPosRankSumfilter" \
  -o "$OUTPUT_DIR"/VariantCalls_fltrdINDELs.vcf

java -cp $GATK37 \
  org.broadinstitute.gatk.tools.CatVariants \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$OUTPUT_DIR"/VariantCalls_fltrdSNPs.vcf \
  -V "$OUTPUT_DIR"/VariantCalls_fltrdINDELs.vcf \
  -out "$OUTPUT_DIR"/VariantCalls_tempcombn.vcf

~/develop/jdk1.8.0_91/jre/bin/java -jar ~/develop/picard/picard.jar \
  SortVcf \
  I="$OUTPUT_DIR"/VariantCalls_tempcombn.vcf \
  O="$OUTPUT_DIR"/VariantCalls_tempsorted.vcf

grep -v ",\*" \
  "$OUTPUT_DIR"/VariantCalls.vcf \
  > "$OUTPUT_DIR"/VariantCalls_temp.vcf
mv "$OUTPUT_DIR"/VariantCalls_temp.vcf "$OUTPUT_DIR"/VariantCalls.vcf
grep -v ",\*" \
  "$OUTPUT_DIR"/VariantCalls_tempsorted.vcf \
  > "$OUTPUT_DIR"/VariantCalls_fltrd.vcf

java -jar $GATK37 \
  -T ValidateVariants \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$OUTPUT_DIR"/VariantCalls.vcf
java -jar $GATK37 \
  -T ValidateVariants \
  -R "$INDEX_DIR"/ConsensusSeqs.fasta \
  -V "$OUTPUT_DIR"/VariantCalls_fltrd.vcf

# Remove temp files
rm "$OUTPUT_DIR"/VariantCalls_SNPsOnly.*
rm "$OUTPUT_DIR"/VariantCalls_fltrdSNPs.*
rm "$OUTPUT_DIR"/VariantCalls_INDELsOnly.*
rm "$OUTPUT_DIR"/VariantCalls_fltrdINDELs.*
rm "$OUTPUT_DIR"/VariantCalls_tempcombn.*
rm "$OUTPUT_DIR"/VariantCalls_tempsorted.*
