#!/usr/bin/env Rscript
library(VariantAnnotation)
library(GenomicRanges)
library(IRanges)

vcf.file <- c("/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/4Founders.deepVar.vcf.gz")
vcf.index.file <- c("/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/4Founders.deepVar.vcf.gz.tbi")

# Read in GM marker annotations
gm_metadata <- read.csv("data/gm_uwisc_v1.csv")
gm_metadata <- gm_metadata[!is.na(gm_metadata$bp_mm10),]
gm_metadata <- gm_metadata[which(gm_metadata$chr != "M"),]
gm_metadata <- data.frame(gm_metadata$chr,
                          gm_metadata$bp_mm10,
                          gm_metadata$bp_mm10,
                          gm_metadata$marker)
colnames(gm_metadata) <- c("chr","start","end","marker")

# Make GRanges object
gm_metadata_GR <- GenomicRanges::GRanges(seqnames = as.character(gm_metadata$chr), 
                                         ranges = IRanges::IRanges(start = gm_metadata$start,
                                                                   end = gm_metadata$end),
                                         marker = gm_metadata$marker,
                                         marker.pos = gm_metadata$start)
supplied_gm_metadata <- data.frame(gm_metadata_GR)
param <- VariantAnnotation::ScanVcfParam(which=gm_metadata_GR)
vcf <- VariantAnnotation::readVcf(vcf.index.file,genome = "mm10",param)

# Pull marker info and ref/alt data
vcf_meta <- rowRanges(vcf)
vcf.meta.2 <- data.frame(vcf_meta)

# Pull genotypes from VCF
genos <- geno(vcf)$GT

# Saving necessary components that can be joined in a tidyverse enviroment
save(vcf.meta.2, 
     genos, 
     supplied_gm_metadata, 
     file = "data/04_pull_4WC_founder_genos_out.RData")

