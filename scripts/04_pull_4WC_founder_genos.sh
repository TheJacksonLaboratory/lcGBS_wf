#!/bin/bash
#SBATCH -J QTL2_VarAnn
#SBATCH --mem 100GB
#SBATCH -t 8:00:00

singularity run docker://sjwidmay/lcgbs_hr:variantannotation scripts/04_pull_4WC_founder_genos.R