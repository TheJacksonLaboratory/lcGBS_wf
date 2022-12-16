#!/bin/bash
#SBATCH -J CC_sam2bam_test
#SBATCH --mem 100GB
#SBATCH -t 24:00:00

fastq_dir=/fastscratch/widmas
cd ${fastq_dir}

singularity run /projects/omics_share/meta/containers/quay.io-biocontainers-samtools-1.14--hb421002_0.img samtools view -S -b test.sam > test.bam
singularity run /projects/omics_share/meta/containers/quay.io-biocontainers-samtools-1.14--hb421002_0.img samtools sort test.bam -o test.sorted.bam
singularity run /projects/omics_share/meta/containers/quay.io-biocontainers-samtools-1.14--hb421002_0.img samtools index test.sorted.bam
