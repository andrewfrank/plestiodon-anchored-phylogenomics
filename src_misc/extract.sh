#!/bin/bash
#SBATCH --job-name=extract
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=andrew.frank@uconn.edu
#SBATCH -o extract_%j.out
#SBATCH -e extract_%j.err

tar -xf /home/CAM/afrank/PlestiodonAnchoredPhylogenomics.tar
tar -xf /home/CAM/afrank/PlestiodonAnchoredPhylogenomics_archive.tar
tar -xf /home/CAM/afrank/WholeGenomeAssemblies_SkiltonianusGilberti.tar
tar -xf /home/CAM/afrank/WholeGenomeAssemblies_SkiltonianusGilberti_archive.tar