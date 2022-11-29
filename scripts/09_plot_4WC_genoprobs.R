#!/usr/bin/env Rscript
library(qtl2)
library(dplyr)
library(tidyr)
library(parallel)
library(qtlcharts)
library(broman)
library(ggbeeswarm)

# Load genotype probabilities and cross data
load("data/4WC_cross.RData")
load("data/4WC_genoprobs.RData")

# A = CAST
# B = POHN
# C = GOR
# D = PWD
# Phase genotypes
ph <- qtl2::guess_phase(cross = X4WC_cross, geno = m)

# Generate Plots
for(i in 1:nrow(X4WC_cross$cross_info)){
  png(file=paste0("plots/plot_onegeno_",rownames(X4WC_cross$cross_info)[i],".png"))
  qtl2::plot_onegeno(geno = ph, 
                     map = map,
                     ind = i, 
                     col = c(qtl2::CCcolors[6],"lightgreen","darkred",qtl2::CCcolors[7])) # add legends here
  dev.off()
}

