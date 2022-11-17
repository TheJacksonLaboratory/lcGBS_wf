#!/usr/bin/env Rscript
library(VariantAnnotation)
library(GenomicRanges)

# Import 4WC vcf if it hasn't been imported already (takes roughly 40 min)
if(!file.exists("data/4WC_vcf.Rdata")){
  # Read in VCF using bioconductor/VariantAnnotation
  X4WC_vcf <- VariantAnnotation::readVcf("/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/4Founders.deepVar.vcf.gz")
  
  # Save this file just in case so it can be loaded
  save(X4WC_vcf, file = "data/4WC_vcf.Rdata")
}

# Load in 4WC vcf
#load("data/4WC_vcf.Rdata")

# Read in GM marker annotations
gm_metadata <- read.csv("data/gm_uwisc_v1.csv")

gm_metadata <- gm_metadata[!is.na(gm_metadata$bp_mm10),]
GenomicRanges::GRanges(seqnames = "chr",
                       start = "bp_mm10",
                       end = "bp_mm10")