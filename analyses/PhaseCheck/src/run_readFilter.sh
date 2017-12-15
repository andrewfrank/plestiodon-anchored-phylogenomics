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

#$ -N samtools
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"
IDNAME="${ID[$SGE_TASK_ID - 1]}"

# Specify directories
DATA_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/bwa_mappings"
BED_DIR="/home/$USER/PlestiodonAnchoredPhylogenomics/PhaseCheck/data"
OUTPUT_DIR="/tempdata3/$USER/PlestiodonAnchoredPhylogenomics/data/filtered_mappings"

# Create output directory if not present
if [ ! -d "$OUTPUT_DIR" ]; then
mkdir -p "$OUTPUT_DIR"
fi

# Initiate an initial filter from the original mapping, saving only reads that
# are properly mapped (both mates mapped to same reference), that are primary
# alignments, and that are not duplicates
#   (-h) keep sam header information
#   (-q 1) keep reads with mapping quality >= 1 (toss unmapped reads)
#   (-f 2) keep reads mapped in proper pairs
#   (-F 256) toss reads not in a primary alignment
samtools view \
  $DATA_DIR/"$DATAFILENAME"_mapping.bam \
  -h \
  -q 1 \
  -f 2 \
  -F 256 \
  -F 1024 \
  -o $OUTPUT_DIR/"$DATAFILENAME"_temp.sam

# Use awk to save headings from the sam file
#   $1 is the first column (used for header info, always kept)
awk \
  '$1 ~ /^@/ {print}' \
  $OUTPUT_DIR/"$DATAFILENAME"_temp.sam \
  > $OUTPUT_DIR/"$DATAFILENAME"_temp.headings

# Use awk to select particular reads based on position and CIGAR string values.
# This goes through all reads line by line, and does filtering per refseq.
# The filter keeps reads that DO NOT have soft clips, keeps reads that are only
# soft clipped at the start of the refseq, and keeps reads that are only soft
# clipped at the end of the refseq.
#   $3 is the refseq column
#   $4 is the POS (position) column
#   $6 is the CIGAR string column
#   ^[^S]+$ selects strings without S (where S indicates soft clipping)
#   ^.+S$ selects strings that might be soft clipped in the beginning and are
#       always soft clipped at the end
#   ^[^S]+S$ selects strings that only end with soft clipping
grep "@SQ" $OUTPUT_DIR/"$DATAFILENAME"_temp.headings | while read -r LINES; do
   REF=$(echo "$LINES" | cut -d $'\t' -f 2 | cut -d ':' -f 2)
   LEN=$(echo "$LINES" | cut -d $'\t' -f 3 | cut -d ':' -f 2)
   awk -v ref="$REF" -v len="$LEN" \
    ' $3 ~ ref &&
    ( $6 ~ /^[^S]+$/ ||
    ( $4 == 1 && $6 !~ /^.+S$/ ) ||
    ( $4 == (len - $6) + 1 && $6 ~ /^[^S]+S$/ ) ) {print}' \
    $OUTPUT_DIR/"$DATAFILENAME"_temp.sam \
    >> $OUTPUT_DIR/"$DATAFILENAME"_fltrdtemp.sam
done

# Go thru samtools view again, and apply the regional excessive coverage filter
# (keep everything that overlaps with the areas indicated in the BED file)
cat \
  $OUTPUT_DIR/"$DATAFILENAME"_temp.headings \
  $OUTPUT_DIR/"$DATAFILENAME"_fltrdtemp.sam |
samtools view \
  - \
  -L "$BED_DIR"/excessive_coverage_filter.bed \
  -U $OUTPUT_DIR/"$DATAFILENAME"_tossedreads.sam \
  -o $OUTPUT_DIR/"$DATAFILENAME"_keptreads.sam \

# Take mates in the tossed reads and see if they have mates in the kept reads,
# then put matched reads back into the kept reads sam
join \
  <(sort -k1,1 $OUTPUT_DIR/"$DATAFILENAME"_keptreads.sam) \
  <(sort -k1,1 $OUTPUT_DIR/"$DATAFILENAME"_tossedreads.sam) \
  -o 0,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15 |
sed "s/ /\t/g" > $OUTPUT_DIR/"$DATAFILENAME"_orphanreads.sam

# Put headings back on this filtered sam file, run samtools sort on the final
# kept reads sam
cat \
  $OUTPUT_DIR/"$DATAFILENAME"_temp.headings \
  $OUTPUT_DIR/"$DATAFILENAME"_keptreads.sam \
  $OUTPUT_DIR/"$DATAFILENAME"_orphanreads.sam \
  > $OUTPUT_DIR/"$DATAFILENAME"_fltrdmapping.sam
samtools sort \
  $OUTPUT_DIR/"$DATAFILENAME"_fltrdmapping.sam \
  -T $OUTPUT_DIR/"$DATAFILENAME"_temp \
  -o $OUTPUT_DIR/"$DATAFILENAME"_fltrdtemp2.bam

# Rerun mark groups on new filtered mapping
~/develop/jdk1.8.0_91/jre/bin/java -jar ~/develop/picard/picard.jar \
  AddOrReplaceReadGroups \
  I="$OUTPUT_DIR"/"$DATAFILENAME"_fltrdtemp2.bam \
  O="$OUTPUT_DIR"/"$DATAFILENAME"_fltrdmapping.bam \
  RGID="$IDNAME" \
  RGLB=lib1 \
  RGPL=illumina \
  RGPU="$IDNAME" \
  RGSM="$DATAFILENAME"

# Run Picard tools to make BAM index on filtered BAM
~/develop/jdk1.8.0_91/jre/bin/java -jar ~/develop/picard/picard.jar \
  BuildBamIndex \
  I=$OUTPUT_DIR/"$DATAFILENAME"_fltrdmapping.bam

# Remove intermediates
rm $OUTPUT_DIR/"$DATAFILENAME"_temp.sam
rm $OUTPUT_DIR/"$DATAFILENAME"_temp.headings
rm $OUTPUT_DIR/"$DATAFILENAME"_fltrdtemp.sam
rm $OUTPUT_DIR/"$DATAFILENAME"_keptreads.sam
rm $OUTPUT_DIR/"$DATAFILENAME"_fltrdmapping.sam
rm $OUTPUT_DIR/"$DATAFILENAME"_fltrdtemp2.bam
