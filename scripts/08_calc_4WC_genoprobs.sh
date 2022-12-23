#!/bin/bash
#SBATCH -J lcGBS_calcgenoprobs
#SBATCH --mem 100GB
#SBATCH -t 8:00:00

singularity run docker://sjwidmay/lcgbs_hr:qtl2_et_al scripts/08_calc_4WC_genoprobs.R
