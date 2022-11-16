#!/usr/bin/env Rscript
library(qtl2)
library(qtl2convert)
library(dplyr)
library(vroom)
library(purrr)

# Inputs are original genotype files from 01_geneseek2qtl2.R
geno_files <- list.files("data/all_genos/", pattern = "DO_4WC_geno")

# Read in covariate file
covar <- read.csv("data/DO4WC_covar.csv")

# Identify mice for each cross
DOmice <- covar %>% dplyr::filter(CrossType == "DO")
DOmice %>% 
  dplyr::select(-CrossType) %>%
  vroom::vroom_write(., path = "data/DO_covar.csv", delim = ",")

X4WCmice <- covar %>% dplyr::filter(CrossType == "4WC")
X4WCmice %>% 
  dplyr::select(-CrossType) %>%
  vroom::vroom_write(., path = "data/4WC_covar.csv", delim = ",")

# Rewrite genotype files specific to samples in each cross; founders are different and will eventually have different founder genotype files
separateCrosses <- function(geno_file){
  g <- vroom::vroom(paste0("data/all_genos/",geno_file), skip = 3, col_names = T)
  
  DOg <- g %>%
    dplyr::select(marker,colnames(g)[which(colnames(g) %in% DOmice$SampleID)])
  newfilename <- paste0("data/DO_genos/",gsub(geno_file, pattern = "4WC_", replacement = ""))
  qtl2convert::write2csv(df = DOg, 
                         filename = newfilename,
                         comment = paste(newfilename, 
                                         "genotypes for chr",
                                         strsplit(strsplit(newfilename, split = "geno")[[1]],split = ".csv")[[2]]), 
                         overwrite = T)

  X4WCg <- g %>%
    dplyr::select(marker,colnames(g)[which(colnames(g) %in% X4WCmice$SampleID)])
  newfilename <- paste0("data/4WC_genos/",gsub(geno_file, pattern = "DO_", replacement = ""))
  qtl2convert::write2csv(df = X4WCg,
                         filename = newfilename,
                         comment = paste(newfilename, 
                                         "genotypes for chr",
                                         strsplit(strsplit(newfilename, split = "geno")[[1]][[3]],split = ".csv")[[1]]), 
                         overwrite = T)

}
purrr::map(geno_files, separateCrosses)

# Write control file for DO mice
chr <- c(1:19, "X")
write_control_file(output_file = "data/DOforqtl2.json",
                   crosstype="do",
                   description="DO_HR",
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("DO_genos/DO_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   covar_file="DO_covar.csv",
                   sex_covar="Sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="Generation", overwrite = T)

# Write control file for 4WC mice
write_control_file(output_file = "data/4WCforqtl2.json",
                   crosstype="genail4",
                   description="4WC_HR",
                   # Update lines below when these consensus genotypes actually exist; reading in this cross throws a warning
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"), 
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("4WC_genos/4WC_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   covar_file="4WC_covar.csv",
                   sex_covar="Sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="Generation", overwrite = T)



