#!/bin/bash
#SBATCH -J CC_fastq_download
#SBATCH --mem 100GB
singularity exec docker://rocker/r-base:latest scripts/download_ena_data.R
