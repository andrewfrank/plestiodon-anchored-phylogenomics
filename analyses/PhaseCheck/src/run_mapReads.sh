#!/bin/bash

## FIXES TO MAKE: INTEGRATE PANDASEQ https://github.com/neufeld/pandaseq, THEN
# SICKLE TRIM, AND THEN RUN BWA-MEM MULTIPLE TIMES, THEN MERGE THE SAM FILES

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

ID=(
H9NDFADXX.2.GGTAGTTG
H9NDFADXX.2.GGTCGGCG
H9NDFADXX.2.GGTTCAGC
H9NDFADXX.2.TCGTTCGA
H9NDFADXX.2.TGACGCGG
H9NDFADXX.2.TGGAGTAC
H9NDFADXX.2.TGGTATGA
H9NDFADXX.2.TTACTTCT
H9NDFADXX.2.TTATAACT
H9NDFADXX.2.TTCGATGA
H9NDFADXX.2.TTGACTAG
HCHLNBCXX.2.TTCGTCGG
HCHLNBCXX.2.TTGGAGCT
HCHLNBCXX.2.TTGGCTCC
)

#$ -N bwa-mem
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"
IDNAME="${ID[$SGE_TASK_ID - 1]}"

# Specify directories
INDEX_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data/bwa-index"
DATA_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/trimmed_reads"
OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/bwa_mappings"

# Create output directory if not present
if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

# Run bwa mem, samtools sort into BAM file, remove original SAM
bwa mem \
  -M \
  "$INDEX_DIR"/ConsensusSeqs \
  $DATA_DIR/"$DATAFILENAME"_R1_trim.fastq \
  $DATA_DIR/"$DATAFILENAME"_R2_trim.fastq \
  > "$OUTPUT_DIR"/"$DATAFILENAME"_tempmapping.sam
samtools sort \
  "$OUTPUT_DIR"/"$DATAFILENAME"_tempmapping.sam \
  -O bam \
  -T "$OUTPUT_DIR"/"$DATAFILENAME"_tempmapping \
  -o "$OUTPUT_DIR"/"$DATAFILENAME"_tempmapping.bam
rm "$OUTPUT_DIR"/"$DATAFILENAME"_tempmapping.sam

# Run Picard Tools, mark duplicates, remove temp mapping intermediate
~/develop/jdk1.8.0_91/jre/bin/java -jar ~/develop/picard/picard.jar \
  MarkDuplicates \
  TAGGING_POLICY=All \
  I="$OUTPUT_DIR"/"$DATAFILENAME"_tempmapping.bam \
  O="$OUTPUT_DIR"/"$DATAFILENAME"_tempmark.bam \
  M="$OUTPUT_DIR"/"$DATAFILENAME"_dupMetrics.txt
rm "$OUTPUT_DIR"/"$DATAFILENAME"_tempmapping.bam

# Run Picard Tools, add read groups, remove temp marked mapping intermediate
~/develop/jdk1.8.0_91/jre/bin/java -jar ~/develop/picard/picard.jar \
  AddOrReplaceReadGroups \
  I="$OUTPUT_DIR"/"$DATAFILENAME"_tempmark.bam \
  O="$OUTPUT_DIR"/"$DATAFILENAME"_mapping.bam \
  RGID="$IDNAME" \
  RGLB=lib1 \
  RGPL=illumina \
  RGPU="$IDNAME" \
  RGSM="$DATAFILENAME"
rm "$OUTPUT_DIR"/"$DATAFILENAME"_tempmark.bam

# Run Picard tools to make BAM index on BAM
~/develop/jdk1.8.0_91/jre/bin/java -jar ~/develop/picard/picard.jar \
  BuildBamIndex \
  I=$OUTPUT_DIR/"$DATAFILENAME"_mapping.bam
