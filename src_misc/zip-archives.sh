#!/bin/bash

#$ -N zip-archives
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -t 1-1:1
#$ -tc 1
#$ -pe smp 1
#$ -S /bin/bash
#$ -o $JOB_NAME_$TASK_ID.out
#$ -e $JOB_NAME_$TASK_ID.err
#$ -cwd

#$ -S /bin/bash

cd /archive/afrank

gzip < PlestiodonAnchoredPhylogenomics.tar > PlestiodonAnchoredPhylogenomics.tar.gz
gzip < PlestiodonAnchoredPhylogenomics_archive.tar > PlestiodonAnchoredPhylogenomics_archive.tar.gz
gzip < WholeGenomeAssemblies_SkiltonianusGilberti.tar > WholeGenomeAssemblies_SkiltonianusGilberti.tar.gz
gzip < WholeGenomeAssemblies_SkiltonianusGilberti_archive.tar > WholeGenomeAssemblies_SkiltonianusGilberti_archive.tar.gz
