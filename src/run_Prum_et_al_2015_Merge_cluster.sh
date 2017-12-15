#!/bin/bash

DATA_FILES=(
# AMF026_Sample_P0077_EJ_I18066
- CAS223944_Sample_P0077_EJ_I7487 #nothing happened?
- EULA_02_Sample_P0077_EJ_I7488 #threw java error
- JQR025_Sample_P0077_EJ_I7489 #ran out of memory
# JQR517_Sample_P0077_EJ_I7490
# JQR552_Sample_P0077_EJ_I7491
- JQR576_Sample_P0077_EJ_I7492 #threw java error
# JQR645_Sample_P0077_EJ_I18068
# JQR792_Sample_P0077_EJ_I7493
- JQR1255_Sample_P0077_EJ_I7495 #threw java error
- JQR1589_Sample_P0077_EJ_I7496 #nothing happened?
- JQR1615_Sample_P0077_EJ_I7497 #ran out of memory
- JQR1870_Sample_P0077_EJ_I18069 #ran out of memory
- SK_02_Sample_P0077_EJ_I7498 #nothing happened?
)

PROJECT_SUBDIR="mtDNA_project"
#$ -N mtDNA_project
#$ -t 1-14
#$ -tc 10

#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

#$ -S /bin/bash
DATAFILENAME="${DATA_FILES[$SGE_TASK_ID - 1]}"

WORKING_DIR="/home/$USER/$PROJECT_SUBDIR"

# Merge.java requires the original fastq file - move into the sample directory
# and unzip the fastq.gz if not already done
cd "$WORKING_DIR/$DATAFILENAME"
for FILE in ./*.fastq.gz
do
  BASE=$(basename $FILE .fastq.gz)  # Basename must have ext written out
  FASTQ="$BASE.fastq"               # Filename for the temporary fastq file
  if [ ! -f $FASTQ ]; then
    gzip -dc $FILE > $FASTQ         # Unzip the fastq
  fi
done

# Get current sample number, then rename the sample directory to the name
# compatible with Merge.java (which is I[SampleNumber])
SAMPLENUMBER=$(echo "$DATAFILENAME" | cut -dI -f 2)
mv "$WORKING_DIR/$DATAFILENAME" \
  "$WORKING_DIR/I$SAMPLENUMBER"

# Merge.java requires you to be in the Code directory to run
cd "$WORKING_DIR/Code"

# Make sure Merge.java is complied, then run
javac Merge.java
java -Xmx2g Merge $SAMPLENUMBER $SAMPLENUMBER 8

# Rename the sample directory back to its original name
mv "$WORKING_DIR/I$SAMPLENUMBER" \
  "$WORKING_DIR/$DATAFILENAME"

# Move back to the sample directory to gzip resultant files and remove generated
# fastq files
cd "$WORKING_DIR/$DATAFILENAME"
gzip "I$SAMPLENUMBER_M.fastq"
gzip "I$SAMPLENUMBER_U1.fastq"
gzip "I$SAMPLENUMBER_U2.fastq"

# Remove any superfluous fastq files from the sample directory
rm *.fastq
