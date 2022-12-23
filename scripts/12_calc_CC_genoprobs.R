#!/usr/bin/env Rscript
library(qtl2)
library(dplyr)
library(ggplot2)
library(parallel)
library(ggbeeswarm)
library(future)
library(furrr)
library(qtl2fst)
pullCrossovers <- function(xo_locations, xo_chrs){
  ind_crossovers <- list()
  for(i in names(xo_locations)){
    if(length(xo_locations[[i]]) == 0){
      next
      ind_crossovers[[i]]
    } else {
      chr_xos <- data.frame(names(xo_locations[i]), data.frame(xo_locations[[i]])) %>%
        `colnames<-`(c("sample","cM")) %>%
        dplyr::mutate(chr = xo_chrs)
      ind_crossovers[[i]] <- chr_xos
    }
  }
  Reduce(rbind, ind_crossovers) %>%
    dplyr::select(sample, chr, cM)
}
gm_meta_build38 <- read.csv("data/gm_uwisc_v1.csv")
gm_meta_build39 <- read.csv("data/gm_uwisc_v2.csv")

#####
# Load in CC cross data
#####
load("data/CC_cross.RData")


if(!file.exists("data/CC_pr_X.fst")){
  #####
  # Insert pseudomarkers
  #####
  map <- qtl2::insert_pseudomarkers(CC_cross$gmap, step = 1)
  save(map, file = "data/CC_map.RData")
  
  #####
  # Calculate genotype probabilities
  #####
  # pr <- qtl2::calc_genoprob(CC_cross, map, error_prob=0.002, cores = parallel::detectCores())
  fpr <- qtl2fst::calc_genoprob_fst(cross = CC_cross, 
                                    map = map, 
                                    error_prob = 0.002, 
                                    cores = parallel::detectCores(), 
                                    fdir = "data/", 
                                    fbase = "CC_pr")
  
  #####
  # Convert geno probs to allele probs
  #####
  # save(pr, file = "data/CC_genoprobs.RData")
  
  #####
  # Find best marginal genotype probability
  #####
  m_char <- qtl2::maxmarg(probs = fpr,
                          cores = parallel::detectCores(),
                          return_char = T,
                          minprob = 0.5)
  save(m_char, file = "data/CC_maxmarg.RData")
  m <- qtl2::maxmarg(probs = pr, 
                     cores = parallel::detectCores(), 
                     minprob = 0.5)
  save(m, file = "data/CC_maxmarg_numeric.RData")
} else {
  load("data/CC_map.RData")
  # load("data/CC_genoprobs.RData")
  load("data/CC_maxmarg.RData")
  load("data/CC_maxmarg_numeric.RData")
}



#####
# Count crossovers
#####

percent_missing <- n_missing(CC_cross, "ind", "prop")*100
nxo <- qtl2::count_xo(m, cores=parallel::detectCores())
nxo_df <- data.frame(nxo) %>%
  `colnames<-`(colnames(nxo)) %>%
  dplyr::mutate(sample = rownames(nxo)) %>%
  dplyr::left_join(., CC_cross$covar %>%
                     dplyr::mutate(sample = rownames(CC_cross$covar)))
long_nxo_df <- nxo_df %>%
  tidyr::pivot_longer(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X),
                      names_to = "chromosome", values_to = "nxos") %>%
  dplyr::mutate(chromosome = as.factor(chromosome))
long_nxo_df$chromosome <- factor(long_nxo_df$chromosome, levels = c(paste(seq(1:19)),"X"))

#####
# Locate crossovers
#####
xolocs <- qtl2::locate_xo(geno = m, 
                          map = map,
                          cores=parallel::detectCores())
# xolocs_tr <- purrr::transpose(xolocs)
# inds <- names(xolocs_tr)
xolocs_chr <- names(xolocs)
future::plan(multisession, workers = 16)
make_chunks <- furrr:::make_chunks
make_chunks(n_x = length(xolocs_chr), n_workers = 16)

# pullCrossovers(xo_locations = xolocs[[3]],
#                xo_chrs = xolocs_chr[[3]])
xolocs_df <- furrr::future_map2(.x = xolocs, 
                                .y = xolocs_chr,
                                .f = pullCrossovers,
                                .options = furrr_options(seed = TRUE)) %>% 
  Reduce(dplyr::bind_rows,.)

xolocs_df_nest <- xolocs_df %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::group_by(chr, sample) %>%
  tidyr::nest()


#####
# Identify diplotypes on both sides of crossover
#####
# pullXoDiplotypes(xo_chrs = xolocs_df_nest$chr[[1]],
#                  xo_dat = xolocs_df_nest$data[[1]],
#                  xo_sample = xolocs_df_nest$sample[[1]])

pullXoDiplotypes <- function(xo_chrs, xo_dat, xo_sample){
  # samples <- unique(xo_dat$sample)
  # sample_xo_diplotype_list <- list()
  cat(paste0(xo_sample, "; (Chromosome ",xo_chrs,")"))
  cMsearch <- data.frame(CC_cross$gmap[xo_chrs]) %>%
    dplyr::mutate(marker = rownames(.)) %>%
    `colnames<-`(c("cM","marker"))
  rownames(cMsearch) <- NULL
  
  dip_list <- apply(xo_dat[,1], 1, function(x){
    if(which.min(abs(cMsearch$cM-x)) < 21){
      search_index <- seq(from = 1, to = 20, by = 1)
    } else if(which.min(abs(cMsearch$cM-x)) > nrow(cMsearch)-21){
      search_index <- rep(which.min(abs(cMsearch$cM-x)),21)-seq(from = 20, to = 0, by = -1)
    } else {
      search_index <- c(rep(which.min(abs(cMsearch$cM-x)),21)-seq(from = 10, to = -10, by = -1))
    }
    cMsearch_df <- cMsearch[search_index,]
    search_pos <- cMsearch_df[!is.na(cMsearch_df$cM),]$cM
    marker <- cMsearch_df[!is.na(cMsearch_df$cM),]$marker
    options(scipen = 99999)
    dips <- matrix(nrow = length(unique(marker)), ncol = 4)
    fpr <- readRDS("data/CC_pr_fstindex.rds")
    
    fpr_sample_chr <- qtl2fst::subset_fst_genoprob(x = fpr, ind = xo_sample, chr = xo_chrs)
    rm(fpr)
    #[ind,chr,mar]
    for(p in unique(marker)){
      
      interval_genoprobs <- qtl2::pull_genoprobpos(genoprobs = fpr_sample_chr, 
                                                   map = CC_cross$gmap,
                                                   chr = xo_chrs,
                                                   marker = p)
      
      
      
      # interval_genoprobs <- qtl2::pull_genoprobpos(genoprobs = fpr,
      #                                              map = CC_cross$gmap,
      #                                              chr = xo_chrs,
      #                                              pos = p)
      probs <- interval_genoprobs[which(interval_genoprobs == max(interval_genoprobs))]
      diplotype <- dimnames(interval_genoprobs)[[2]][which(interval_genoprobs == max(interval_genoprobs))]
    
      # write diplotype
      dips[which(unique(marker) == p),1] <- diplotype
      dips[which(unique(marker) == p),2] <- search_pos[which(unique(marker) == p)]
      dips[which(unique(marker) == p),3] <- probs
      dips[which(unique(marker) == p),4] <- p
    }
    dips <- data.frame(dips)
    colnames(dips) <- c("diplotype","cM","prob","marker")
    for(i in 2:nrow(dips)){
      if(dips[i-1,]$diplotype == dips[i,]$diplotype){
        next
      } else {
        dips <- dips[c(i-1,i),] %>%
          dplyr::mutate(loc = c("proximal","distal"))
        break
      }
    }
    return(dips)
  })
  
  dip_list[[which(unique(samp$cM) == cM)]] <- dips
    sample_xo_diplotype_list[[which(samples == s)]] <- Reduce(dplyr::bind_rows, dip_list) %>%
      dplyr::mutate(sample = s) %>%
      dplyr::select(sample, everything())
  return(Reduce(dplyr::bind_rows, sample_xo_diplotype_list) %>%
           dplyr::mutate(chr = xo_chrs))
}

future::plan(multisession, workers = 16)
make_chunks <- furrr:::make_chunks
make_chunks(n_x = length(xolocs_df_nest$data), n_workers = 16)
xodiplotypes <- furrr::future_map2(.x = xolocs_df_nest$chr, 
                                   .y = xolocs_df_nest$data,
                                   .f = pullXoDiplotypes, 
                                   .options = furrr_options(seed = TRUE))
all_xo_diplotypes <- list()
for(i in 1:length(xodiplotypes)){
  all_xo_diplotypes[[i]] <- xodiplotypes[[i]] %>%
    dplyr::filter(!is.na(loc))
}
xos_sample_nested <- Reduce(rbind,all_xo_diplotypes) %>%
  dplyr::group_by(sample) %>%
  tidyr::nest()

#####
# Write crossover summary table
#####
auto_xo_count_df <- purrr::map2_dfr(.x = xos_sample_nested$sample, 
                                    .y = xos_sample_nested$data, 
                                    .f = function(smp, data){
                                      data %>%
                                        dplyr::filter(loc == "proximal",
                                                      chr != "20") %>%
                                        dplyr::distinct() %>%
                                        dplyr::count() %>%
                                        dplyr::rename(auto_xos = n) %>% 
                                        dplyr::mutate(sample = smp) %>%
                                        dplyr::select(sample, auto_xos)}) %>%
  dplyr::arrange(sample)
write.csv(auto_xo_count_df, file = "data/CC_autosomal_crossovers.csv", row.names = F, quote = F)

#####
# Plot the number of crossovers
#####
nxo_ggplot <- CC_cross$covar %>%
  dplyr::mutate(sample = rownames(.)) %>%
  dplyr::left_join(.,auto_xo_count_df) %>%
  dplyr::select(sample, auto_xos, Generation) %>%
  ggplot(., mapping = aes(x = reorder(sample, auto_xos), y = auto_xos, colour = Generation)) +
  theme_bw() + 
  geom_point() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Sample",
       y = "Number of Autosomal Crossovers")
ggsave(nxo_ggplot, filename = "plots/CC_nxos.png",
       height = 7, width = 7)


#####
# Write individual sample crossover diplotypes
#####
dir.create(path = "data/CC_crossover_diplotypes/", showWarnings = F)
purrr::map2(.x = xos_sample_nested$data, 
            .y = xos_sample_nested$sample, 
            .f = function(x,y){
              outfile <- x %>%
                dplyr::mutate(sample = y) %>%
                dplyr::left_join(., gm_meta_build39, 
                                 by = c("chr", "marker")) %>% 
                dplyr::distinct(sample, chr, cM, marker, loc, prob, bp_mm10, bp_grcm39) %>%
                dplyr::arrange(chr)
              write.csv(outfile,
                        file = paste0("data/CC_crossover_diplotypes/",y,"_crossover_diplotypes.csv"), 
                        row.names = F, quote = F)
})

