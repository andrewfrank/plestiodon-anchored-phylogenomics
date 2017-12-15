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

PROJECT_SUBDIR="mapping_project"
#$ -N ph_pipeline
#$ -t 1-14
#$ -tc 20
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

WORKING_DIR="/home/$USER/$PROJECT_SUBDIR/$DATAFILENAME"
if [ ! -d "$WORKING_DIR" ]; then
mkdir -p "$WORKING_DIR"
fi
cd "$WORKING_DIR"

SCRATCH_DIR="/scratch/afrank/$PROJECT_SUBDIR/$DATAFILENAME"
if [ ! -d "$SCRATCH_DIR" ]; then
mkdir -p "$SCRATCH_DIR"
fi
cd "$SCRATCH_DIR"

# Establish variables for file names
## Data prep filenames
RUN_STATUS="$WORKING_DIR/status.txt"

FQR1_CAT=$DATAFILENAME.R1cat.fq.gz
FQR2_CAT=$DATAFILENAME.R2cat.fq.gz

## Trim Galore output filenames
FQR1_TRIMMED=$DATAFILENAME.R1cat_val_1.fq
FQR2_TRIMMED=$DATAFILENAME.R2cat_val_2.fq

# Concatenate read files
echo "Concatenating..." > $RUN_STATUS
cat $WORKING_DIR/*R1*.fastq.gz > ./$FQR1_CAT
cat $WORKING_DIR/*R2*.fastq.gz > ./$FQR2_CAT

echo "File sizes:" >> $RUN_STATUS
echo `wc -c $FQR1_CAT` >> $RUN_STATUS
echo `wc -c $FQR2_CAT` >> $RUN_STATUS

# Run trim galore to clean reads, output to fastq
echo "Trimming..." >> $RUN_STATUS
trim_galore \
  --paired $FQR1_CAT $FQR2_CAT \
  --phred33 \
  --illumina \
  --dont_gzip \
  --output_dir $SCRATCH_DIR
# --output_dir $WORKING

echo "File sizes:" >> $RUN_STATUS
echo `wc -c $FQR1_TRIMMED` >> $RUN_STATUS
echo `wc -c $FQR2_TRIMMED` >> $RUN_STATUS

for FILE in /home/$USER/$PROJECT_SUBDIR/fasta_files/T176_L*.fasta
do

  # Establish variables for file names
  ## Bowtie2 filenames
  BASENAME=$(basename $FILE .fasta)
  INDEX_BASENAME=$DATAFILENAME.$BASENAME
  INDEX_FILTER=$INDEX_BASENAME.txt
  INDEX_FA=$INDEX_BASENAME.fa
  SAM=$INDEX_BASENAME.sam
  BAM=$INDEX_BASENAME.bam

  ## samtools filenames
  TMP_BAM_SORT=$INDEX_BASENAME.tmp
  BAM_SORT=$INDEX_BASENAME.sorted.bam
  VCF=$INDEX_BASENAME.vcf
  HC=$INDEX_BASENAME.HapCompass

  # Get the first sequence from the matching individual from the alignment fasta
  grep "^>$DATAFILENAME" $FILE |
    cut -d'>' -f2 |
    head -n 1 > $INDEX_FILTER

  # Filter this loci's fasta file according current individual
  echo "Filtering $FILE multi-fasta by individual of interest..." >> $RUN_STATUS
  seqtk subseq $FILE $INDEX_FILTER > $INDEX_FA
  echo "File sizes:" >> $RUN_STATUS
  echo `wc -c $INDEX_FA` >> $RUN_STATUS

  # Build Bowtie2 index for this loci's data
  echo "Creating $INDEX_BASENAME bowtie2 index..." >> $RUN_STATUS
  bowtie2-build -f --seed 12345 $INDEX_FA $INDEX_BASENAME

  # Run bowtie2
  echo "Running bowtie2 for $INDEX_BASENAME..." >> $RUN_STATUS
  bowtie2 \
    -x "$INDEX_BASENAME" \
    -q \
    --phred33 \
    --sensitive-local \
    -I 250 \
    -X 800 \
    --dovetail \
    -1 $FQR1_TRIMMED \
    -2 $FQR2_TRIMMED \
    -S $SAM
  echo "File sizes:" >> $RUN_STATUS
  echo `wc -c $SAM` >> $RUN_STATUS

  # Convert SAM to BAM, and keep ONLY uniquely mapped reads
  # from: https://www.biostars.org/p/56246/
  echo "Filtering $SAM..." >> $RUN_STATUS
  samtools view -bq 1 $SAM -o $BAM
  rm $SAM
  echo "File sizes:" >> $RUN_STATUS
  echo `wc -c $BAM` >> $RUN_STATUS

  # Sort and index bam files
  echo "Sorting $SAM..." >> $RUN_STATUS
  samtools sort -T $TMP_BAM_SORT -o $BAM_SORT $BAM
  echo "File sizes:" >> $RUN_STATUS
  echo `wc -c $BAM_SORT` >> $RUN_STATUS

  # Modify fasta reference files - change all ambi. bases to a stnd base
  ## samtools mpileup doesn't work with IACUC ambiguous bases in the reference,
  ## due to this issue: https://github.com/samtools/hts-specs/issues/54
  echo "Removing ambiguous bases from $INDEX_FA..." >> $RUN_STATUS
  sed -i '$ s/Y/C/g;$ s/R/A/g;$ s/W/A/g;$ s/S/G/g;$ s/K/T/g;$ s/M/C/g;$ s/D/A/g;$ s/V/A/g;$ s/H/A/g;$ s/B/C/g;' $INDEX_FA

  # Run samtools mpileup and bcftools call
  echo "Running samtools mpileup / bcftools call on $BAM_SORT..." >> $RUN_STATUS
  samtools mpileup \
    --skip-indels \
    -uf $INDEX_FA \
    $BAM_SORT |
  bcftools call \
    -mv \
    --output $VCF
  echo "File sizes:" >> $RUN_STATUS
  echo `wc -c $VCF` >> $RUN_STATUS

  if [[ $BASENAME == T176_L24 ]] || [[ $BASENAME == T176_L33 ]] || [[ $BASENAME == T176_L336 ]]
  then
    # Run HapCompass allowing for tetraploidy, since these loci are paralogs
    echo "Running HapCompass on $BAM_SORT and $VCF..." >> $RUN_STATUS
    ~/develop/jdk1.8.0_91/bin/java -jar \
      ~/develop/hapcompass_v0.8.2/hapcompass.jar \
      --ploidy 4 \
      -b $BAM_SORT \
      -v $VCF \
      -o $HC/$BASENAME
  else
    # Run HapCompass assuming diploidy
    echo "Running HapCompass on $BAM_SORT and $VCF..." >> $RUN_STATUS
    ~/develop/jdk1.8.0_91/bin/java -jar \
      ~/develop/hapcompass_v0.8.2/hapcompass.jar \
      -b $BAM_SORT \
      -v $VCF \
      -o $HC/$BASENAME
  fi

  # Prep files for return to home directory
  echo "Moving files back to home directory..." >> $RUN_STATUS
  cp ./$BAM $WORKING_DIR/$BAM_SORT
  cp ./$INDEX_FA $WORKING_DIR/$INDEX_FA
  cp ./$VCF $WORKING_DIR/$VCF
  cp -R ./$HC $WORKING_DIR/$HC

done

## Clean up scratch directory
echo "Removing scratch directory..." >> $RUN_STATUS
rm -r $SCRATCH_DIR/*
