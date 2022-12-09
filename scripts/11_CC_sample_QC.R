#!/usr/bin/env Rscript
library(qtl2)
library(dplyr)

CC_cross <- qtl2::read_cross2("data/cc_qtl2_genail.json")
CC_cross <- qtl2::drop_nullmarkers(CC_cross)

# Reordering genotypes so that most common allele in founders is first

for(chr in seq_along(CC_cross$founder_geno)) {
  fg <- CC_cross$founder_geno[[chr]]
  g <- CC_cross$geno[[chr]]
  f1 <- colSums(fg==1)/colSums(fg != 0)
  
  fg[fg==0] <- NA
  g[g==0] <- NA
  
  fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
  g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
  
  fg[is.na(fg)] <- 0
  g[is.na(g)] <- 0
  
  CC_cross$founder_geno[[chr]] <- fg
  CC_cross$geno[[chr]] <- g
}

# percent missing genotypes - CC
percent_missing <- qtl2::n_missing(CC_cross, "ind", "prop")*100
missing_genos_df <- data.frame(names(percent_missing), percent_missing) %>%
  `colnames<-`(c("sample","percent_missing"))
missing_genos_plot <- ggplot(data = missing_genos_df, mapping = aes(x = sample, 
                                                                    y = percent_missing)) + 
  theme_bw() + 
  geom_point(shape = 21) + 
  ylim(c(0,100)) +
  labs(title = "DO Missing Genotypes") + 
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(missing_genos_plot, filename = "plots/CC_missing_genos.png", width = 6, height = 6)
save(CC_cross, file = "data/CC_cross.RData")
