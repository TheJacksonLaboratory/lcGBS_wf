#!/bin/bash
#SBATCH -J CC_bowtie2_test
#SBATCH --mem 100GB
#SBATCH -t 24:00:00

fastq_dir=/fastscratch/widmas
ref_genome_base=/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bowtie2/Mus_musculus.GRCm38.dna.primary_assembly.fa
cd ${fastq_dir}
singularity run docker://quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0 bowtie2 -x ${ref_genome_base} -1 test_forward_paired.fq.gz -2 test_reverse_paired.fq.gz -S test.sam
