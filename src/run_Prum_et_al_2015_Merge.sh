#!/bin/bash

WORKING_DIR="/home/$USER/$PROJECT_SUBDIR/Code"
cd "$WORKING_DIR"

SAMPLENUMBER=$(echo "$DATAFILENAME" | cut -dI -f 2)
echo $SAMPLENUMBER

javac Merge.java
java -Xmx2g Merge $SAMPLENUMBER $SAMPLENUMBER 8
