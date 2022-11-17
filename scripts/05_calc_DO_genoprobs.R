#!/usr/bin/env Rscript
library(qtl2)
library(dplyr)
library(ggplot2)
library(parallel)

# Load in DO cross data
load("data/DO_cross.RData")

# Insert pseudomarkers
map <- qtl2::insert_pseudomarkers(DO_cross$gmap, step = 1)

# Calculate genotype probabilities
pr <- calc_genoprob(DO_cross, map, error_prob=0.002, cores = parallel::detectCores())

# Convert geno probs to allele probs
save(pr, file = "data/DO_genoprobs.RData")