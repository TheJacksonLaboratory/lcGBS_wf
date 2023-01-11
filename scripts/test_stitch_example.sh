#!/bin/bash
#SBATCH -J STITCH_TEST
#SBATCH --mem 3000GB
#SBATCH -t 1:30:00
#SBATCH -p high_mem

singularity run docker://sjwidmay/lcgbs_hr:stitch /projects/compsci/vmp/USERS/widmas/lcGBS_wf/scripts/test_STITCH.R
