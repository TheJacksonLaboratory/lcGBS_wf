#!/usr/bin/env Rscript
library(STITCH)
library(parallel)

# Input files
mouse_datadir <- file.path("test/wgs/mouse/CC_subsampled_fastqs")

# Path to original bamlist
mouse_bamlist <- file.path(mouse_datadir, "bams/STITCH_bamlist.txt")
temp_bamlist <- read.table(mouse_bamlist)
# mouse_genfile <- paste0(mouse_datadir, "gen.txt")

# Create sample name file in lieu of reheadering .bams
# Avoiding reheadering .bams because initial testing and validation may 
# use simulated or altered data
sample_names <- gsub(unlist(lapply(strsplit(temp_bamlist[,1], 
                                            split = "/"), 
                                   function(x) x[[7]])),
                     pattern = "_dedup.bam", replacement = "")
write.table(x = sample_names, 
            file = "sample_names.txt",
            quote = F, 
            row.names = F, 
            col.names = F,
            sep = "\t")
sample_names_file <- "sample_names.txt"

# Path to pos file
# pos file is a headerless txt file with four columns:
# 1) chr_ where _ is a number
# 2) position
# 3) ref allele
# 4) alt allele

mouse_posfile <- "/fastscratch/STITCH_outputDir/work/b7/9f3f7c5a990b6ddc39449bcbb73efd/STITCH_9_pos.txt"
temp_posfile <- read.table(mouse_posfile)

# Remove duplicate positions
temp_posfile <- temp_posfile[!duplicated(temp_posfile$V2),]

# Write a temporary pos file to be called by STITCH
write.table(x = temp_posfile, 
            file = "temp_pos.txt",
            quote = F, 
            row.names = F, 
            col.names = F,
            sep = "\t")


# Specify this file for the function below
dup_removed_posfile <- "temp_pos.txt"

# Set number of founders
# This should be a param in nextflow
mouse_K <- 8

# Number of generations
# This can be a param in nextflow, but maybe not necessary
# mouse_nGen <- 100

# Chromosome to be analyzed
mouse_chr <- "chr9"

STITCH::STITCH(tempdir = tempdir(),
               chr = mouse_chr,
               bamlist = mouse_bamlist,
               posfile = dup_removed_posfile,
               outputdir = paste0(getwd(), "/"),
               sampleNames_file = sample_names_file,
               K = mouse_K, 
               nGen = 20)