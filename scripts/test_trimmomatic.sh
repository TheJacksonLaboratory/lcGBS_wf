#!/bin/bash
#SBATCH -J CC_trim_test
#SBATCH --mem 100GB
#SBATCH -t 8:00:00

fastq_dir=/fastscratch/widmas
cd ${fastq_dir}
singularity run docker://quay.io/biocontainers/trimmomatic:0.35--6 trimmomatic PE GES15-06696-CC001-Unc_CCGTCCCG_HNLK2CCXX_L001_001.R1.fastq.gz GES15-06696-CC001-Unc_CCGTCCCG_HNLK2CCXX_L001_001.R2.fastq.gz test_forward_paired.fq.gz test_forward_unpaired.fq.gz test_reverse_paired.fq.gz test_reverse_unpaired.fq.gz LEADING:3 TRAILING:3 MINLEN:36
