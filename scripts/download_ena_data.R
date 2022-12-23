#!/usr/bin/env Rscript
options(timeout = max(3600, getOption("timeout"))) # change timeout time in seconds, default is 60. 3600 is 1 hour.
cc_data = read.csv("data/cc_genome_download/CCStrains.csv")
cc_ena_files = read.table("data/cc_genome_download/filereport_read_run_PRJEB14673_tsv.txt",
                          sep = "\t", header = TRUE)

cc_ena_fastqs = cc_ena_files[grep("\\.fastq\\.gz$",cc_ena_files$submitted_ftp),]
cc_ena_bams = cc_ena_files[grep("\\.fastq\\.gz$",cc_ena_files$submitted_ftp,invert=TRUE),]

# Turning FASTQ files index into useful output
cc_ena_fastqs$fastq_1_ftp = gsub("^(.*);(.*)$","\\1",cc_ena_fastqs$fastq_ftp)
cc_ena_fastqs$fastq_2_ftp = gsub("^(.*);(.*)$","\\2",cc_ena_fastqs$fastq_ftp)
cc_ena_fastqs$submitted_1_ftp = gsub("^(.*);(.*)$","\\1",cc_ena_fastqs$submitted_ftp)
cc_ena_fastqs$submitted_2_ftp = gsub("^(.*);(.*)$","\\2",cc_ena_fastqs$submitted_ftp)
cc_ena_fastqs$stripped_1_ftp = gsub("ftp.sra.ebi.ac.uk/vol1/run/ERR\\d{3}/ERR\\d+/(.*\\.fastq.gz)$","\\1",cc_ena_fastqs$submitted_1_ftp)
cc_ena_fastqs$stripped_2_ftp = gsub("ftp.sra.ebi.ac.uk/vol1/run/ERR\\d{3}/ERR\\d+/(.*\\.fastq.gz)$","\\1",cc_ena_fastqs$submitted_2_ftp)
cc_ena_fastqs$cc_strain = gsub("^GES15-\\d{5}-(.*)-(.*)_[ATCG]{8}_.*$","\\1/\\2",cc_ena_fastqs$stripped_1_ftp)

# Turning BAM files index into useful output
cc_ena_bams$strain = gsub("^ftp.sra.ebi.ac.uk/vol1/run/ERR\\d{3}/ERR\\d+/(.*)_(.*)_[MF]\\d+_[NYGCU]{3,4}\\.sorted\\.bam$","\\1/\\2",cc_ena_bams$submitted_ftp)

# # Downloading files
for(i in 1:nrow(cc_ena_files)){
  fastqs = unlist(strsplit(cc_ena_files[i,"submitted_ftp"], ";"))
  fastqs_stripped = gsub("ftp.sra.ebi.ac.uk/vol1/run/ERR191/ERR\\d+/(.*\\.fastq.gz)$","\\1",fastqs)
  for(j in 1:length(fastqs)){
    download.file(fastqs[j], paste0("/fastscratch/widmas/",fastqs_stripped[j]))
  }
}
