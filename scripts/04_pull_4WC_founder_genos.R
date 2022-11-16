#!/usr/bin/env Rscript
library(VariantAnnotation)
library(GenomicRanges)

X4WC_vcf <- VariantAnnotation::readVcf("/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/4Founders.deepVar.vcf.gz")
save(X4WC_vcf, file = "data/4WC_vcf.Rdata")
