#!/usr/bin/env Rscript
library(qtl2)
library(qtl2fst)
library(purrr)
library(tidyr)
library(lsa)
library(ggplot2)
library(cowplot)
library(data.table)
library(tictoc)
setwd("/projects/compsci/vmp/USERS/widmas/stitch-nf/")


## GATK data
# project directory where original fastqs were
fastqDir <- "data/DO_ddRADseq_NovaSeq/"

# make plotting directory
plotDir <- paste0(fastqDir,"cossim_plots")
dir.create(path = plotDir, showWarnings = F)

# # function to find path to sample qtl2 files
# find_gbs_genoprobs <- function(sample, fastqDir){
#   sampleDir <- list.files(paste0(fastqDir,"qtl2files/"), pattern = sample)
#   return(paste0(fastqDir,"qtl2files/",sampleDir))
# }

# find genoprobs
# genoprobDir <- find_gbs_genoprobs(sample, fastqDir)
genoprobDir <- paste0(fastqDir,"qtl2files")
genoprobs <- list.files(genoprobDir, pattern = "36_state")

# load genoprobs
load(file = paste(genoprobDir,genoprobs,sep = "/"))

# load allele probs
alleleprobs <- list.files(genoprobDir, pattern = "8_state")

# load genoprobs
load(file = paste(genoprobDir,alleleprobs,sep = "/"))

# view cross
cross

# markers missing genotypes in x samples
options(scipen = 9999)
hist(qtl2::n_missing(cross, by = "marker"), 
     main = "number of markers missing genotypes for x samples", 
     xlab = "x samples")

# samples missing genotypes in x markers
options(scipen = 9999)
hist(qtl2::n_missing(cross, by = "individual"),
     main = "number of samples missing genotypes for x markers", 
     xlab = "x markers")

# make GATK viterbi
GATK_m <- qtl2::maxmarg(probs = fpr,
                        cores = (parallel::detectCores()/2),
                        minprob = 0.01)

# Phase genotypes
GATK_ph <- qtl2::guess_phase(cross = cross,
                             geno = GATK_m)

### GigaMUGA data
# load original geno probs & cross data
load("../lcGBS_wf/data/DO_cross.RData")
load("../lcGBS_wf/data/DO_genoprobs.RData")
apr_DO_cross <- genoprob_to_alleleprob(pr, cores = (parallel::detectCores()/2), quiet = F)
load("../lcGBS_wf/data/DO_maxmarg_numeric.RData")

# Phase genotypes
ph <- qtl2::guess_phase(cross = DO_cross,
                             geno = m)

# create a map between sample names and GM sample names
sample_key <- list()
for(k in 1:nrow(DO_cross$covar)){
  gbs_sample <- rownames(cross$covar)[grep(pattern = rownames(DO_cross$covar[k,]), 
                                           x = rownames(cross$covar))]
  if(length(gbs_sample) == 0){
    next
  }
  GM_sample <- rownames(DO_cross$covar[k,])
  sample_key[[k]] <- data.frame(gbs_sample, GM_sample)
}
sample_key_df <- purrr::discard(.x = sample_key, .p = is.null) %>%
  Reduce(rbind,.)

# make chromosome plots
# for(sam in sample(x = seq(1:nrow(sample_key_df)),size = 5)){
#   # plot GM probs
#   qtl2::plot_onegeno(geno = ph,
#                      map = DO_cross$pmap,
#                      ind = sample_key_df[sam,]$GM_sample,
#                      col = c(qtl2::CCcolors),
#                      main = paste0("GigaMUGA - ",sample_key_df[sam,]$GM_sample)) # add legends here
#   
#   # plot GATK probs
#   qtl2::plot_onegeno(geno = GATK_ph,
#                      map = cross$pmap,
#                      ind = sample_key_df[sam,]$gbs_sample,
#                      main = paste0("GATK - ",sample_key_df[sam,]$GM_sample))
# }

# Calculate cosine similarity between GM allele probs and
# GATK allele probs
sample_cossim_data <- list()
for(sam in 1:nrow(sample_key_df)){
  print(sample_key_df[sam,]$GM_sample)
  cossim_all_chrs <- list()
  cossim_founders_all_chrs <- list()
  for(c in 1:19){
  
  print(c)
  # find GM marker positions
  gm_pos <- unique(DO_cross$pmap[[c]])
  names(gm_pos) <- NULL

  # find GATK marker positions
  gatk_pos <- unique(cross$pmap[[c]])
  names(gatk_pos) <- NULL

  # find concordant marker sets
  GATK_index <- list()
  GM_index <- list()
  for(i in 1:length(gm_pos)){
  foo <- gatk_pos[which(abs(gm_pos[i]-gatk_pos) == min(abs(gm_pos[i]-gatk_pos)))][1]
  
  if(i != 1){
    if(foo == max(unlist(GATK_index))){
      next
    } else {
      GATK_index[[i]] <- foo
      GM_index[[i]] <- gm_pos[i]
    }
  } else {
    GATK_index[[i]] <- foo
    GM_index[[i]] <- gm_pos[i]
  }
}

  gatk_matching <- gatk_pos[gatk_pos %in% unlist(GATK_index)]
  gm_matching <- gm_pos[gm_pos %in% unlist(GM_index)]

  matching_pmap_gatk <- cross$pmap[[c]][cross$pmap[[c]] %in% gatk_matching]
  matching_pmap_gm <- unique(DO_cross$pmap[[c]])[unique(DO_cross$pmap[[c]]) %in% gm_matching]
  foo <- DO_cross$pmap[[c]][DO_cross$pmap[[c]] %in% matching_pmap_gm]

  gm_common_sites <- names(foo[!duplicated(foo)])
  gatk_common_sites <- names(matching_pmap_gatk)

  options(scipen = 99)
  cossim_list_markers <- list()
  # sampled_samples <- 1:nrow(sample_key_df)
  
  tictoc::tic()
  for(s in 1:length(gm_common_sites)){
  # for(s in sample(1:length(gm_common_sites), size = 100, replace = F)){
      cossim <- lsa::cosine(x = apr_DO_cross[[c]][sample_key_df[sam,]$GM_sample,,gm_common_sites[s]],
                            y = apr[[c]][sample_key_df[sam,]$gbs_sample,,gatk_common_sites[s]])
  
      char_allele <- sort(apr_DO_cross[[c]][sample_key_df[sam,]$GM_sample,,gm_common_sites[s]], decreasing = T)[1:2]
      soft_diplotype <- char_allele[sort(names(char_allele))]
      soft_diplotype <- data.frame(t(data.frame(c(soft_diplotype, names(soft_diplotype)))))
      rownames(soft_diplotype) <- NULL
      colnames(soft_diplotype) <- c("pA1","pA2","A1","A2")
      cossim_samples <- cbind(data.frame(matching_pmap_gatk[s], foo[!duplicated(foo)][s], cossim, sample_key_df[sam,]$GM_sample), soft_diplotype)
      cossim_list_markers[[s]] <- cossim_samples
  }
  cossim_df_markers <- Reduce(rbind, cossim_list_markers)
  tictoc::toc()
  
  # dplyr::mutate(CHROM = c)
                # pA1 = as.numeric(pA1),
                # pA2 = as.numeric(pA2))
  colnames(cossim_df_markers) <- c("GATK_pos","GigaMUGA_pos", "CosSim","sample","pA1","pA2","A1","A2")
  cossim_df_markers$CHROM <- c
  rownames(cossim_df_markers) <- NULL
  
  cossim_df_markers[which(cossim_df_markers$pA1 > 0.85),]$A2 <- cossim_df_markers[which(cossim_df_markers$pA1 > 0.85),]$A1
  cossim_df_markers[which(cossim_df_markers$pA2 > 0.85),]$A1 <- cossim_df_markers[which(cossim_df_markers$pA2 > 0.85),]$A2
  cossim_df_markers <- cossim_df_markers %>%
    dplyr::mutate(A1 = as.factor(A1),
                A2 = as.factor(A2))

  color_model <- names(qtl2::CCcolors)
  names(color_model) <- LETTERS[1:8]

  levels(cossim_df_markers$A1) <- names(qtl2::CCcolors)[names(color_model) %in% levels(cossim_df_markers$A1)]
  levels(cossim_df_markers$A2) <- names(qtl2::CCcolors)[names(color_model) %in% levels(cossim_df_markers$A2)]
  cossim_all_chrs[[c]] <- cossim_df_markers
  
  A1 <- ggplot() +
    theme_bw() +
    geom_point(data = cossim_df_markers,
             mapping = aes(x = GATK_pos,
                           y = CosSim,
                           colour = A1),
             alpha = 0.7,
             size = 1) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none") +
    scale_colour_manual(values = qtl2::CCcolors) +
    facet_wrap(.~sample, ncol = 4) +
    ggtitle("A1 Call")+
    ylim(c(0,1))

  A2 <- ggplot() +
    theme_bw() +
    geom_point(data = cossim_df_markers,
               mapping = aes(x = GATK_pos,
                             y = CosSim,
                             colour = A2),
               alpha = 0.7,
               size = 1) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none") +
    scale_colour_manual(values = qtl2::CCcolors) +
    facet_wrap(.~sample, ncol = 4) +
    ggtitle("A2 Call")+
    ylim(c(0,1))

  cossim_list_founders <- list()
  for(f in 1:8){
  cossim <- lsa::cosine(x = apr_DO_cross[[c]][sample_key_df[sam,]$GM_sample,f,gm_common_sites],
                            y = apr[[c]][sample_key_df[sam,]$gbs_sample,f,gatk_common_sites])
  founders_samples <- data.frame(LETTERS[f], cossim,sample_key_df[sam,]$GM_sample)  
  cossim_list_founders[[f]] <- founders_samples 
  }
  cossim_df_founders <- Reduce(rbind, cossim_list_founders) %>%
    dplyr::mutate(CHROM = c)
  colnames(cossim_df_founders) <- c("founder", "CosSim", "sample","CHROM")
  rownames(cossim_df_founders) <- NULL
  cossim_founders_all_chrs[[c]] <- cossim_df_founders

  byfounder <- ggplot() +
    theme_bw() +
    geom_jitter(data = cossim_df_founders,
                mapping = aes(x = founder,
                              y = CosSim), width = 0.1) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ggtitle(label = paste0("Chromosome ", c)) +
    ylim(c(0,1))

  cossim_scan <- cowplot::plot_grid(A1, A2, ncol = 1)
  all_plots <- cowplot::plot_grid(cossim_scan, byfounder, ncol = 2, rel_widths = c(1,0.6))

  ggsave(all_plots,
       filename = paste0(plotDir,"/",sample_key_df[sam,]$GM_sample,"_chr",c,".png"),
       width = 7, height = 5)
}
  cossim_all_df <- Reduce(rbind, cossim_all_chrs)
  sample_cossim_data[[sam]] <- cossim_all_df
  
  GATK_histogram <- ggplot() + 
    theme_bw() + 
    geom_histogram(data = cossim_all_df, mapping = aes(x = CosSim), bins = 50) + 
    facet_wrap(.~CHROM) + 
    theme(panel.grid.minor = element_blank(), 
          axis.text.x = element_text(size = 9)) + 
    labs(x = "Cosine Similarity Between GM and GATK Allele Probability Vectors") + 
    ggtitle(sample_key_df[sam,]$GM_sample, subtitle = "GATK")
  ggsave(GATK_histogram, 
         filename = paste0(plotDir,"/",sample_key_df[sam,]$GM_sample,"_cossim.png"), width = 8, height = 8)
}

sample_cossim_data_df <- Reduce(rbind, sample_cossim_data)
save(sample_cossim_data_df, file = paste0(plotDir,"/GATKCosSim.RData"))
