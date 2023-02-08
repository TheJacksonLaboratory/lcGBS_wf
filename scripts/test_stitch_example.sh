#!/bin/bash
#SBATCH -J STITCH_TEST
#SBATCH --mem 10GB
#SBATCH -t 1:30:00
#SBATCH -p compute

singularity run docker://sjwidmay/stitch_nf:latest /projects/compsci/vmp/USERS/widmas/stitch-nf/bin/stitch/run_stitch.R
