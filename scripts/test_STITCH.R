#!/usr/bin/env Rscript
library(STITCH)
library(parallel)

STITCH::STITCH(tempdir = tempdir(), 
               chr = "chr19", 
               bamlist = "bamlist.txt", 
               posfile = "pos.txt", 
               genfile = "gen.txt", 
               outputdir = paste0(getwd(), "/"), 
               K = 4, 
               nGen = 100, 
               nCores = parallel::detectCores())