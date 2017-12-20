#!/bin/bash

#$ -N BBCtoXanadu
#$ -M andrew.frank@uconn.edu
#$ -m besa
#$ -cwd

cd /tempdata3/afrank

tar -cvf WholeGenomeAssemblies_SkiltonianusGilberti.tar WholeGenomeAssemblies_SkiltonianusGilberti
md5sum WholeGenomeAssemblies_SkiltonianusGilberti.tar
scp ./WholeGenomeAssemblies_SkiltonianusGilberti.tar afrank@xanadu-submit-ext.cam.uchc.edu:/home/CAM/afrank
