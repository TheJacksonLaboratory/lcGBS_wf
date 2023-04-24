library(data.table)
library(qtl2)
library(dplyr)
library(plotly)
library(qtl2fst)
setwd("/projects/compsci/vmp/USERS/widmas/stitch-nf/")
samples <- data.table::fread("data/CC_data/qtl2files/geno10.csv", skip = 3, nrows = 1)
samples <- colnames(samples)[-1]
cross_info <- data.table::fread("data/CC_data/qtl2files/test_crossinfo.csv", nrows = 10) %>%
  dplyr::mutate(id_2 = samples) %>%
  dplyr::select(-id) %>%
  dplyr::select(id_2, everything()) %>%
  dplyr::rename(id = id_2)
write.csv(cross_info, file = "data/CC_data/qtl2files/crossinfo.csv", row.names = F, quote = F)


# Write new control file for CC mice
chr <- unlist(lapply(unique(list.files("data/CC_data/qtl2files/", 
                                pattern = "foundergeno")), 
              function(x) gsub(strsplit(x, ".csv")[[1]][[1]], 
                               pattern = "foundergeno", 
                               replacement = "")))

qtl2::write_control_file(output_file = "data/CC_data/qtl2files/cc_qtl2_genail.json",
                         crosstype="genail8",
                         description="CC_GENAIL",
                         founder_geno_file=paste0("foundergeno", chr, ".csv"),
                         founder_geno_transposed=TRUE,
                         gmap_file=paste0("gmap", chr, ".csv"),
                         pmap_file=paste0("pmap", chr, ".csv"),
                         geno_file=paste0("geno", chr, ".csv"),
                         geno_transposed=TRUE,
                         geno_codes=list(A=1, H=2, B=3),
                         xchr="X",
                         # sex_covar = "Sex",
                         # sex_codes=list(female = "female", male = "male"), 
                         # covar_file="cc_qtl2_genail/test_covar.csv", 
                         crossinfo_file = "crossinfo.csv",
                         overwrite = T)

# Load in the cross object
stitch_CC <- qtl2::read_cross2("data/CC_data/qtl2files/cc_qtl2_genail.json")

# Drop null markers
stitch_CC <- qtl2::drop_nullmarkers(stitch_CC)

for(chr in seq_along(stitch_CC$founder_geno)) {
  fg <- stitch_CC$founder_geno[[chr]]
  g <- stitch_CC$geno[[chr]]
  f1 <- colSums(fg==1)/colSums(fg != 0)
  
  fg[fg==0] <- NA
  g[g==0] <- NA
  
  fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
  g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
  
  fg[is.na(fg)] <- 0
  g[is.na(g)] <- 0
  
  stitch_CC$founder_geno[[chr]] <- fg
  stitch_CC$geno[[chr]] <- g
}

# Calculate the percent of missing genotypes per sample
percent_missing <- qtl2::n_missing(stitch_CC, "ind", "prop")*100
missing_genos_df <- data.frame(names(percent_missing), percent_missing) %>%
  `colnames<-`(c("sample","percent_missing"))


# Plot missing genotypes per sample
missing_genos_plot <- ggplot(data = missing_genos_df, mapping = aes(x = reorder(sample, percent_missing),
                                                                    y = percent_missing)) +
  theme_bw() +
  geom_point(shape = 21) +
  labs(title = "DO Missing Genotypes") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plotly::ggplotly(missing_genos_plot)

# Identify duplicates??
cg <- qtl2::compare_geno(stitch_CC, cores=2)
qtl2::plot_compare_geno(x = cg, rug = T, main = "CC subsampled reads")

# Calculate genotype probs
dir.create("data/CC_data/qtl2files/genoprobs", showWarnings = F)
fpr <- suppressWarnings(qtl2fst::calc_genoprob_fst(cross = stitch_CC,
                                                   map = stitch_CC$gmap,
                                                   fbase = "pr",
                                                   fdir = "data/CC_data/qtl2files/genoprobs",
                                                   error_prob=0.002,
                                                   overwrite=TRUE,
                                                   quiet = F,
                                                   cores = 16))
save(fpr, file = "data/CC_data/qtl2files/genoprobs_36_state.RData")

# calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = fpr)
dir.create("data/CC_data/qtl2files/alleleprobs", showWarnings = F)
fapr <- qtl2fst::fst_genoprob(genoprob = apr,
                              fbase = "apr", 
                              fdir = "data/CC_data/qtl2files/alleleprobs",
                              overwrite=TRUE,
                              quiet = F)
save(fapr, file = "data/CC_data/qtl2files/genoprobs_8_state.RData")

# find maxmarg
m <- qtl2::maxmarg(probs = fpr, 
                   cores = 16, 
                   minprob = 0.8,
                   quiet = F)

# Phase genotypes
ph <- qtl2::guess_phase(cross = stitch_CC, geno = m)

# Generate Plots

for(i in 1:nrow(stitch_CC$cross_info)){
  # png(file=paste0("data/CC_data/qtl2files/plots/haplo_",rownames(stitch_CC$cross_info)[i],".png"))
  qtl2::plot_onegeno(geno = ph, 
                     map = stitch_CC$pmap,
                     ind = i, 
                     col = c(qtl2::CCcolors), main = rownames(stitch_CC$cross_info)[i]) # add legends here
  # dev.off()
}


# read CC001 diplotypes
CC001_dips <- data.table::fread("/projects/omics_share/mouse/GRCm38/supporting_files/cc/genoprobs/CC001_Uncb38V01.csv")
CC068_dips <- data.table::fread("/projects/omics_share/mouse/GRCm38/supporting_files/cc/genoprobs/CC068_TauUncb38V01.csv")

mchar <- qtl2::maxmarg(probs = fpr, 
                   cores = 16, 
                   minprob = 0.002, 
                   quiet = F, return_char = T)
# begin comparing GM diplotypes to stitch diplotypes
CC001_dips %>%
  dplyr::filter(chromosome == "11")

data.frame(t(mchar$`11`[1:3,]))




