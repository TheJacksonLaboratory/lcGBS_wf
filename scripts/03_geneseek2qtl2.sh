#!/bin/bash
#SBATCH -J GS2QTL2
#SBATCH --mem 100GB
#SBATCH -t 8:00:00

# Performing Marker QC on GigaMUGA samples
singularity run docker://sjwidmay/lcgbs_hr:qtl2_et_al scripts/03_geneseek2qtl2.R
