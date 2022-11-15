#!/usr/bin/env Rscript
library(qtl2)
chr <- c(1:19, "X")
write_control_file(output_file = "data/DO4WC_forqtl2.json",
                   crosstype="do",
                   description="DO4WC_HR",
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("DO_4WC_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   covar_file="data/DO4WC_covar.csv",
                   sex_covar="Sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="Generation", overwrite = T)
