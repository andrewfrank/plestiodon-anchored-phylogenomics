#!/bin/bash

SAMPLES=(I7487 I7488 I7489 I7490 I7491 I7492 I7493 I7495 I7496 I7497 I7498
I18066 I18068 I18069)

SCRIPT_PATH="/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/analyses/PhaseCheck/src/phasedalign-to-vcf.py"
WORKING_DIR="/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/data"
OUTPUT_DIR="/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/data/lemmons_vcfs"

# Runs my python script to create the VCF, then removes non-variant bases, then
# makes IUPAC ambiguous bases in ALT into symbolic alleles for use in downstream
# software (GATK), then converts leading "n" to a dash, then removes base
# positions containing a dash or a stand-alone N
for SAMPLE in "${SAMPLES[@]}"; do
  REF="$WORKING_DIR"/consensus_seqs/concat-consensus_seqs_noAmbi.fasta
  AL1="$WORKING_DIR"/cleaned_seqs_byIndividual/"$SAMPLE"_seq1.fasta
  AL2="$WORKING_DIR"/cleaned_seqs_byIndividual/"$SAMPLE"_seq2.fasta
  OUT="$OUTPUT_DIR"/"$SAMPLE"_lemmons-tmp.vcf

  python $SCRIPT_PATH $SAMPLE $REF $AL1 $AL2 $OUT

  fgrep -v -e "0/0" "$OUTPUT_DIR"/"$SAMPLE"_lemmons-tmp.vcf |
  sed -e '6,$s/M/<M>/' \
      -e '6,$s/R/<R>/' \
      -e '6,$s/W/<W>/' \
      -e '6,$s/S/<S>/' \
      -e '6,$s/Y/<Y>/' \
      -e '6,$s/K/<K>/' \
      -e '6,$s/V/<V>/' \
      -e '6,$s/H/<H>/' \
      -e '6,$s/D/<D>/' \
      -e '6,$s/B/<B>/' |
  sed -e '6,$s/n/-/3' |
  fgrep -v -e "-" -e "$(printf '\t')N$(printf '\t')" > "$OUTPUT_DIR"/"$SAMPLE"_lemmons.vcf
  rm $OUT
done
