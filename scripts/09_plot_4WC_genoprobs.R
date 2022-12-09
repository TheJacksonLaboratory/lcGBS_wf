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
load("data/4WC_maxmarg_numeric.RData")
load("data/4WC_map.RData")

# A = CAST
# B = POHN
# C = GOR
# D = PWD
# Phase genotypes
ph <- qtl2::guess_phase(cross = X4WC_cross, geno = m)

X4WCcolors <- c(qtl2::CCcolors[6],
  "#5ADBFF",
  "#153B50",
  qtl2::CCcolors[7])
names(X4WCcolors)[2:3] <- c("POHN","GOR")

# Generate Plots
for(i in 1:nrow(X4WC_cross$cross_info)){
  png(file=paste0("plots/plot_onegeno_",rownames(X4WC_cross$cross_info)[i],".png"))
  qtl2::plot_onegeno(geno = ph, 
                     map = map,
                     ind = i, 
                     col = X4WCcolors) # add legends here
  legend(17, 85, 
         legend=c("CAST", "POHN", "GOR", "PWD"),
         fill=X4WCcolors, cex=0.8)
  dev.off()
}