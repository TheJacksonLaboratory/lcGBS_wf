#!/usr/bin/env Rscript
library(qtl2)
library(dplyr)
library(ggplot2)
library(qtlcharts)
library(broman)
library(fst)
DO_cross <- qtl2::read_cross2("data/DOforqtl2.json")
DO_cross <- qtl2::drop_nullmarkers(DO_cross)

X4WC_cross <- qtl2::read_cross2("data/4WCforqtl2.json")
X4WC_cross <- qtl2::drop_nullmarkers(X4WC_cross)

# Reordering genotypes so that most common allele in founders is first
for(chr in seq_along(DO_cross$founder_geno)) {
  fg <- DO_cross$founder_geno[[chr]]
  g <- DO_cross$geno[[chr]]
  f1 <- colSums(fg==1)/colSums(fg != 0)
  
  fg[fg==0] <- NA
  g[g==0] <- NA
  
  fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
  g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
  
  fg[is.na(fg)] <- 0
  g[is.na(g)] <- 0
  
  DO_cross$founder_geno[[chr]] <- fg
  DO_cross$geno[[chr]] <- g
}
for(chr in seq_along(X4WC_cross$founder_geno)) {
  fg <- X4WC_cross$founder_geno[[chr]]
  g <- X4WC_cross$geno[[chr]]
  f1 <- colSums(fg==1)/colSums(fg != 0)
  
  fg[fg==0] <- NA
  g[g==0] <- NA
  
  fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
  g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
  
  fg[is.na(fg)] <- 0
  g[is.na(g)] <- 0
  
  X4WC_cross$founder_geno[[chr]] <- fg
  X4WC_cross$geno[[chr]] <- g
}

# percent missing genotypes - DO
percent_missing <- qtl2::n_missing(DO_cross, "ind", "prop")*100
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
ggsave(missing_genos_plot, filename = "plots/DO_missing_genos.png", width = 6, height = 6)

# percent missing genotypes - 4WC
percent_missing_4WC <- qtl2::n_missing(X4WC_cross, "ind", "prop")*100
missing_genos_df <- data.frame(names(percent_missing_4WC), percent_missing_4WC) %>%
  `colnames<-`(c("sample","percent_missing"))
missing_genos_plot <- ggplot(data = missing_genos_df, mapping = aes(x = sample, 
                                                                    y = percent_missing_4WC)) + 
  theme_bw() + 
  geom_point(shape = 21) + 
  ylim(c(0,100)) + 
  labs(title = "4WC Missing Genotypes") + 
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
missing_genos_plot
ggsave(missing_genos_plot, filename = "plots/4WC_missing_genos.png", width = 6, height = 6)

## DO Sex checks
## Reading in all probe intensities
int <- fst::read.fst("data/DO4WC_intensities.fst")
int <- int[seq(1, nrow(int), by=2),-(1:2)] + int[-seq(1, nrow(int), by=2),-(1:2)]
DOint <- int[,which(colnames(int) %in% qtl2::ind_ids(DO_cross))]
X4WCint <- int[,which(colnames(int) %in% qtl2::ind_ids(X4WC_cross))]

# Interactive plot of DO array intensities per sample
n <- names(sort(percent_missing, decreasing=TRUE))
DO_iboxplot <- iboxplot(log10(t(DOint[,n])+1), orderByMedian=FALSE, chartOpts=list(ylab="log10(SNP intensity + 1)"))
save(DO_iboxplot, file = "plots/DO_array_interactiveboxplot.RData")
# Interactive plot of 4WC array intensities per sample
n <- names(sort(percent_missing_4WC, decreasing=TRUE))
X4WC_iboxplot <- iboxplot(log10(t(X4WCint[,n])+1), orderByMedian=FALSE, chartOpts=list(ylab="log10(SNP intensity + 1)"))
save(X4WC_iboxplot, file = "plots/4WC_array_interactiveboxplot.RData")

## Reading in sex chromosome intensities
xint <- qtl2::read_csv_numer(filename = "data/all_genos/DO_4WC_chrXint.csv")
yint <- qtl2::read_csv_numer(filename = "data/all_genos/DO_4WC_chrYint.csv")
DO_covar <- read.csv("data/DO_covar.csv")
DOxint <- xint[,colnames(xint)[which(colnames(xint) %in% DO_covar$SampleID)]]
DOyint <- yint[,colnames(yint)[which(colnames(yint) %in% DO_covar$SampleID)]]
sex <- substr(colnames(DOxint), 1, 1)

## Testing for uninformative markers and filtering those out
x_pval <- apply(DOxint, 1, function(a) t.test(a ~ sex)$p.value)
y_pval <- apply(DOyint, 1, function(a) t.test(a ~ sex)$p.value)
DOxint_ave <- colMeans(DOxint[x_pval < 0.05/length(x_pval),], na.rm=TRUE)
DOyint_ave <- colMeans(DOyint[y_pval < 0.05/length(y_pval),], na.rm=TRUE)

## Plotting sex chromosome intensities to verify sexes
xyints <- data.frame(DOxint_ave, DOyint_ave) %>%
  dplyr::mutate(sample = rownames(.)) %>%
  dplyr::left_join(., DO_covar %>% dplyr::rename(sample = SampleID))
rownames(xyints) <- NULL
labels <- paste0(names(DOxint_ave), " (", round(percent_missing), "%)")
DOsexcheck_plot <- ggplot(data = xyints, mapping = aes(x = DOxint_ave, 
                                    y = DOyint_ave, 
                                    fill = Sex, 
                                    label = labels)) +
  theme_bw() + 
  geom_point(shape = 21, size = 4) + 
  scale_fill_manual(values = c("green","purple")) + 
  labs(x = "Average X chr intensity",
       y = "Average Y chr intensity")
plotly::ggplotly(DOsexcheck_plot, tooltip = "label")


## 4WC Sex Checks
## Reading in sex chromosome intensities
X4WC_covar <- read.csv("data/4WC_covar.csv")
X4WCxint <- xint[,colnames(xint)[which(colnames(xint) %in% X4WC_covar$SampleID)]]
X4WCyint <- yint[,colnames(yint)[which(colnames(yint) %in% X4WC_covar$SampleID)]]
X4WCmetadata <- data.frame(colnames(X4WCxint)) %>%
  `colnames<-`(c("sample")) %>%
  dplyr::left_join(., X4WC_covar %>%
                     dplyr::rename(sample = SampleID))
sex <- X4WCmetadata$Sex

## Testing for uninformative markers and filtering those out
x_pval <- apply(X4WCxint, 1, function(a) t.test(a ~ sex)$p.value)
y_pval <- apply(X4WCyint, 1, function(a) t.test(a ~ sex)$p.value)
X4WCxint_ave <- colMeans(X4WCxint[x_pval < 0.05/length(x_pval),], na.rm=TRUE)
X4WCyint_ave <- colMeans(X4WCyint[y_pval < 0.05/length(y_pval),], na.rm=TRUE)

## Plotting sex chromosome intensities to verify sexes
xyints <- data.frame(X4WCxint_ave, X4WCyint_ave) %>%
  dplyr::mutate(sample = rownames(.)) %>%
  dplyr::left_join(., X4WC_covar %>% dplyr::rename(sample = SampleID))
rownames(xyints) <- NULL
labels <- paste0(names(X4WCxint_ave), " (", round(percent_missing_4WC), "%)")
X4WCsexcheck_plot <- ggplot(data = xyints, mapping = aes(x = X4WCxint_ave, 
                                                       y = X4WCyint_ave, 
                                                       fill = Sex, 
                                                       label = labels)) +
  theme_bw() + 
  geom_point(shape = 21, size = 4) + 
  scale_fill_manual(values = c("blue","orange")) + 
  labs(x = "Average X chr intensity",
       y = "Average Y chr intensity")
plotly::ggplotly(X4WCsexcheck_plot, tooltip = "label")


## Sample Duplicates
cg <- compare_geno(DO_cross, cores=0)
summary(cg)
cg <- compare_geno(X4WC_cross, cores=0)
summary(cg)

save(DO_cross, file = "data/DO_cross.RData")
save(X4WC_cross, file = "data/4WC_cross.RData")


## Genotype Frequencies
## DO_cross
g <- do.call("cbind", DO_cross$geno[1:19])
fg <- do.call("cbind", DO_cross$founder_geno[1:19])
g <- g[,colSums(fg==0)==0]
fg <- fg[,colSums(fg==0)==0]
fgn <- colSums(fg==3)
gf_ind <- vector("list", 4)
for(i in 1:4) {
  gf_ind[[i]] <- t(apply(g[,fgn==i], 1, function(a) table(factor(a, 1:3))/sum(a != 0)))
}

png(file=paste0("plots/DOcross_genotype_frequencies.png"))
par(mfrow=c(2,2), mar=c(0.6, 0.6, 2.6, 0.6))
for(i in 1:4) {
  triplot(c("AA", "AB", "BB"), main=paste0("MAF = ", i, "/8"))
  tripoints(gf_ind[[i]], pch=21, bg="lightblue")
  tripoints(c((1-i/8)^2, 2*i/8*(1-i/8), (i/8)^2), pch=21, bg="violetred")
  
  if(i>=3) { # label mouse with lowest het
    wh <- which(gf_ind[[i]][,2] == min(gf_ind[[i]][,2]))
    tritext(gf_ind[[i]][wh,,drop=FALSE] + c(0.02, -0.02, 0),
            names(wh), adj=c(0, 1))
  }
  
  # label other mice
  if(i==1) {
    lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.3]
  }
  else if(i==2) {
    lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.48]
  }
  else if(i==3) {
    lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.51]
  }
  else if(i==4) {
    lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.6]
  }
  
  for(ind in lab) {
    if(grepl("^F", ind) && i != 3) {
      tritext(gf_ind[[i]][ind,,drop=FALSE] + c(-0.01, 0, +0.01), ind, adj=c(1,0.5))
    } else {
      tritext(gf_ind[[i]][ind,,drop=FALSE] + c(0.01, 0, -0.01), ind, adj=c(0,0.5))
    }
  }
}
dev.off()


# ## 4WC ## NOT GENERALIZABLE TO THINGS OTHER THAN DO (probably because different MAF distribution)
# g <- do.call("cbind", X4WC_cross$geno[1:19])
# fg <- do.call("cbind", X4WC_cross$founder_geno[1:19])
# g <- g[,colSums(fg==0)==0]
# fg <- fg[,colSums(fg==0)==0]
# fgn <- colSums(fg==3)
# gf_ind <- vector("list", 4)
# for(i in 1:4) {
#   gf_ind[[i]] <- t(apply(g[,fgn==i], 1, function(a) table(factor(a, 1:3))/sum(a != 0)))
# }
# 
# png(file=paste0("plots/4WC_genotype_frequencies.png"))
# par(mfrow=c(2,2), mar=c(0.6, 0.6, 2.6, 0.6))
# for(i in 1:4) {
#   triplot(c("AA", "AB", "BB"), main=paste0("MAF = ", i, "/8"))
#   tripoints(gf_ind[[i]], pch=21, bg="lightblue")
#   tripoints(c((1-i/8)^2, 2*i/8*(1-i/8), (i/8)^2), pch=21, bg="violetred")
#   
#   if(i>=3) { # label mouse with lowest het
#     wh <- which(gf_ind[[i]][,2] == min(gf_ind[[i]][,2]))
#     tritext(gf_ind[[i]][wh,,drop=FALSE] + c(0.02, -0.02, 0),
#             names(wh), adj=c(0, 1))
#   }
#   
#   # label other mice
#   if(i==1) {
#     lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.3]
#   }
#   else if(i==2) {
#     lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.48]
#   }
#   else if(i==3) {
#     lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.51]
#   }
#   else if(i==4) {
#     lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.6]
#   }
#   
#   for(ind in lab) {
#     if(grepl("^F", ind) && i != 3) {
#       tritext(gf_ind[[i]][ind,,drop=FALSE] + c(-0.01, 0, +0.01), ind, adj=c(1,0.5))
#     } else {
#       tritext(gf_ind[[i]][ind,,drop=FALSE] + c(0.01, 0, -0.01), ind, adj=c(0,0.5))
#     }
#   }
# }
# dev.off()

