#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=STITCH_CosSim
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 5:00:00
#SBATCH --mem=250G
#SBATCH --ntasks=1

containerDir=/projects/compsci/vmp/USERS/widmas/stitch-nf/singularity_cache
singularity exec ${containerDir}/lcgbs_hr_qtl2_et_al.sif Rscript bin/stitch/stitch_allele_probs.R $1