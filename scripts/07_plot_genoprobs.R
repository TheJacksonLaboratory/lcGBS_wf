#!/usr/bin/env Rscript
library(qtl2)
library(dplyr)
library(tidyr)
library(parallel)
library(qtlcharts)
library(broman)
library(ggbeeswarm)

# Load genotype probabilities and cross data
load("data/DO_cross.RData")
load("data/DO_genoprobs.RData")

# Insert pseudomarkers in genetic map
map <- insert_pseudomarkers(DO_cross$gmap, step=1)

# Find best marginal genotype probability
m_char <- qtl2::maxmarg(probs = pr,
                        cores = parallel::detectCores(),
                        return_char = T,
                        minprob = 0.5)
save(m_char, file = "data/DO_maxmarg.RData")

m <- qtl2::maxmarg(probs = pr, 
                   cores = parallel::detectCores(), 
                   minprob = 0.5)

# Number of crossovers
percent_missing <- n_missing(DO_cross, "ind", "prop")*100
nxo <- qtl2::count_xo(m, cores=parallel::detectCores())
nxo_df <- data.frame(nxo) %>%
  `colnames<-`(colnames(nxo)) %>%
  dplyr::mutate(sample = rownames(nxo)) %>%
  dplyr::left_join(., DO_cross$covar %>%
                     dplyr::mutate(sample = rownames(DO_cross$covar)))
long_nxo_df <- nxo_df %>%
  tidyr::pivot_longer(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X), 
                      names_to = "chromosome", values_to = "nxos") %>%
  dplyr::mutate(chromosome = as.factor(chromosome))
long_nxo_df$chromosome <- factor(long_nxo_df$chromosome, levels = c(paste(seq(1:19)),"X"))

# Plot number of crossovers
nxo_ggplot <- ggplot(long_nxo_df, mapping = aes(x = chromosome, y = nxos, fill = Sex)) + 
  theme_bw() + 
  geom_point(shape = 21, alpha = 0.7, position = position_quasirandom()) + 
  theme(panel.grid = element_blank()) + 
  labs(x = "Chromosome",
       y = "Number of Crossovers")
ggsave(nxo_ggplot, filename = "plots/DOcross_nxos.png",
       height = 7, width = 7)


# Phase genotypes
ph <- qtl2::guess_phase(cross = DO_cross, geno = m)

# Generate Plots
for(i in 1:nrow(DO_cross$cross_info)){
  png(file=paste0("plots/plot_onegeno_",rownames(DO_cross$cross_info)[i],".png"))
  qtl2::plot_onegeno(geno = ph, 
                     map = map,
                     ind = i, 
                     col = c(qtl2::CCcolors)) # add legends here
  dev.off()
}

