library(qtl2)
library(dplyr)
library(tidyr)
library(parallel)

# Load genotype probabilities and cross data
load("data/DO_cross.RData")
load("data/DO_genoprobs.RData")

# Insert pseudomarkers in genetic map
map <- insert_pseudomarkers(DO_cross$gmap, step=1)

# Find best marginal genotype probability
m_char <- qtl2::maxmarg(probs = pr,
                        cores = parallel::detectCores(),
                        return_char = T,
                        minprob = 0.00001)
save(m_char, file = "data/DO_maxmarg.RData")

m <- qtl2::maxmarg(probs = pr, 
                   cores = parallel::detectCores(), 
                   minprob = 0.00001)

# Phase genotypes
ph <- qtl2::guess_phase(cross = DO_cross, geno = m)

# Generate Plots
for(i in 1:nrow(DO_cross$cross_info)){
  png(file=paste0("plots/plot_onegeno_",rownames(DO_cross$cross_info)[i],".png"))
  qtl2::plot_onegeno(geno = ph, 
                     map = map,
                     ind = i, 
                     col = c(qtl2::CCcolors))
  dev.off()
}

