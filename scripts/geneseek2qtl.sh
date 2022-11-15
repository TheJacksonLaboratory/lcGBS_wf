#!/bin/bash
#SBATCH -J GS2QTL2
#SBATCH --mem 100GB
#SBATCH -t 8:00:00

# Performing Marker QC on GigaMUGA samples
singularity run docker://sjwidmay/do4wc_hr:latest scripts/geneseek2qtl2.R
