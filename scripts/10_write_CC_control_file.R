#!/usr/bin/env Rscript
library(qtl2)

# Write new control file for CC mice
qtl2::write_control_file(output_file = "data/cc_qtl2_genail.json",
                   crosstype="genail8",
                   description="CC_GENAIL",
                   founder_geno_file=paste0("cc_qtl2_genail/test_foundergeno.csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("cc_qtl2_genail/test_gmap.csv"),
                   pmap_file=paste0("cc_qtl2_genail/test_pmap.csv"),
                   geno_file=paste0("cc_qtl2_genail/test_geno.csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   sex_covar = "Sex",
                   sex_codes=list(female = "female", male = "male"), 
                   covar_file="cc_qtl2_genail/test_covar.csv", 
                   crossinfo_file = "cc_qtl2_genail/test_crossinfo.csv",
                   overwrite = T)