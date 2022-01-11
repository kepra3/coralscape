# Title: The effect of microniches on snp data using ordinations
# Author: Katharine Prata
# Date created: 10/01/21

# Packages
library(vcfR)
library(adegenet)
library(ggplot2)
library(vegan)

# Functions
organise_snp_matrix <- function(genind) {
  # turn into snp matrix and remove NULL loci - replace with most common allele 
  # might need to do this by population or species instead
  snp.mat <- genind$tab
  cat(paste0("Missing data before imputations is: ", sum(is.na(snp.mat))))
  cat("\nDimensions are ")
  cat(paste0(c("Individuals: ", "snps: "), dim(snp.mat)))
  cat("\nReplacing NULL snps with most common allele")
  gen.imp <- apply(snp.mat, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
  cat(paste0("\nMissing data after imputations is: ", sum(is.na(gen.imp))))
return(gen.imp)
}

combine_metadata <- function(gen.subset.imp, struc.complex, clusters) {
  # match
  struc.complex.subset <- subset(struc.complex, X %in% rownames(gen.subset.imp))
  clusters.subset <- subset(clusters, Individual %in% rownames(gen.subset.imp))
  # sort
  sorted.clusters.subset <- clusters.subset[order(clusters.subset$Individual),]
  # check if identiical
  cat("Checking gen imp and clusters: ", identical(rownames(gen.subset.imp),sorted.clusters.subset[,1]))
  cat("\nChecking gen imp and struc complex: ", identical(rownames(gen.subset.imp),struc.complex.subset[,1]))
  cat("\nChecking struc complex and clusters: ", identical(sorted.clusters.subset[,1], struc.complex.subset[,1]))
  # Put metadata together
  meta.info <- sorted.clusters.subset
  meta.info[,20:24] <- struc.complex.subset[,2:6]
  return(meta.info)
}

overhang_barplot <- function(meta.info, factor) {
  overhang <- table(meta.info$overhang, factor) %>% matrix(nrow = 2, ncol = length(levels(factor)))
  colnames(overhang) <- levels(factor)
  rownames(overhang) <- levels(meta.info$overhang)
  barplot(overhang, main = "Presence of overhang", xlab = "Taxa",
          col = c("Red", "Green"), ylim = c(0,100))
  legend("topright", c("No", "Yes"), fill = c("Red", "Green"))
}

height_prop_plot <- function(meta.info, factor) {
  
  levels <- levels(factor)
  sample.sizes <- NA
  
  for (i in 1:length(levels)) {
    sample.sizes[i] <- length(factor[factor == levels[i]])
  }
  
  label = paste0(levels, "\n n = ", sample.sizes)
  
  p <- ggplot(meta.info, aes(x = factor, y = prop)) + 
    geom_boxplot(notch = TRUE) + geom_point(size = 3, shape = 1) + 
    scale_x_discrete(labels = label) + theme_classic() + 
    ggtitle("Height proportion")
  return(p)
}

theta_boxplot <- function(df, factor) {
  
  levels <- levels(factor)
  sample.sizes <- NA
  
  for (i in 1:length(levels)) {
    sample.sizes[i] <- length(factor[factor == levels[i]])
  }
  
  label = paste0(levels, "\n n = ", sample.sizes)
  
  p <- ggplot(df, aes(x = factor, y = theta)) +
    geom_boxplot(notch = TRUE) +
    geom_point(size = 3, shape = 1) +
    ggtitle("Angles of attachment") +
    scale_x_discrete(labels = label) +
    theme_classic(base_size = 12)
  return(p)
}

# Import datasets
struc.complex <- read.csv("~/git/coralscape_open3d/results/sample_metadata.csv")
aga.genind <- vcfR2genind(read.vcfR("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/all-aga_1d_nc_20.vcf"))
aga.clusters <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/all-aga_1d_nc_20_6.csv")

# Organise data
gen.imp <- organise_snp_matrix(aga.genind) # WARNING IMPUTATION METHOD PROBABLY INCORRECT
gen.subset.imp <- subset(gen.imp, rownames(gen.imp) %in% struc.complex$X)
meta.info <- combine_metadata(gen.subset.imp, struc.complex, aga.clusters)

# Clean working space
rm(gen.imp)
rm(aga.genind)
rm(struc.complex)
rm(aga.clusters)

# Data structure
str(meta.info)
meta.info[, c(2:6, 13, 24)] <- data.frame(lapply(meta.info[, c(2:6, 13, 24)], as.factor))
levels(meta.info$Clusters) <- c("AA1", "AL1", "AA2", "AL2") # have to change if using different datasets
str(meta.info)

# Environmental data ####
# Overhang
overhang_barplot(meta.info, meta.info$Species)  # Species
overhang_barplot(meta.info, meta.info$Clusters)  # Cluster

# Theta
meta.info$theta <- NA
for (row in 1:length(meta.info$xz)) {
  if (meta.info[row, 21] > meta.info[row, 22]) {
    meta.info[row, 25] <- meta.info[row, 21]
  } else {
    meta.info[row, 25] <- meta.info[row, 22]
  }
}
theta_boxplot(meta.info, meta.info$Species)
theta_boxplot(meta.info, meta.info$Clusters)
# maybe colony rugosity would be good to include here!

# relative depth proportion species
height_prop_plot(meta.info, meta.info$Species)
height_prop_plot(meta.info, meta.info$Clusters)

#### RDA ####
full.rda <- rda(gen.subset.imp ~ theta * overhang * prop, data = meta.info, scale = T)
RsquareAdj(full.rda) # 0.09%
summary(eigenvals(full.rda, model = "constrained"))
# RDA1 explains 74% of the variations, RDA2 explains 8%. RDA3 4.8%, RDA4 4.1%
summary(eigenvals(full.rda, model = "unconstrained"))
# PC1 explains 0.4% of the variation, PC2 0.4%, PC3 0.4%, PC4 0.4%
summary(full.rda)

#Partitioning of correlations:
#  Inertia Proportion
#Total           21764     1.0000
#Constrained      3867     0.1777
#Unconstrained   17897     0.8223

step(full.rda)
# removing each facotr the AIC decreases and removing clusters the AIC increases
# best model only with clusters
# prop overhang useful
# get rid of all interaction terms

# Significance
signif.full.c <- anova.cca(full.rda)
signif.full.c 
# Model is signicant 0.004

# Reduced rda ####
reduced.rda <- rda(gen.subset.imp ~ theta + overhang + prop, data = meta.info, scale = T)
# Include Site as a interactive term with all values and depth as continuous
RsquareAdj(reduced.rda) # 0.09%
summary(eigenvals(reduced.rda, model = "constrained"))
# RDA1 explains 90% of the variations, RDA2 explains 5% RDA3 4.3%
summary(eigenvals(reduced.rda, model = "unconstrained"))
# PC1 explains 42% of the variation, PC2 4%, PC3 2%, PC4 2%
summary(reduced.rda)

#Partitioning of correlations:
#  Inertia Proportion
#Total           21764     1.0000
#Constrained      2835     0.1303
#Unconstrained   18929     0.8697

step(reduced.rda)
# suggests to remove theta but will keep for now as it only slightly improves <2AIC

# Significance
signif.reduced.rda <- anova.cca(reduced.rda)
signif.reduced.rda
# Model is signicant 0.001
signif.axis.reduced.rda <- anova.cca(reduced.rda, by = "axis")
# Error: vector memory exhausted (limit reached?)
# signif.axis.reduced.rda
# Model is signicant 0.001
signif.term.reduced.rda <- anova.cca(reduced.rda, by = "term")
signif.term.reduced.rda

# Arguments within rda function
# permutation = how(within = Within(type = "factor", mirror = TRUE),
# plots = Plots(strata = Depth, type = "factor"), blocks = Loc))
# Linear model relative
# lme(Y ~ X, random = ~1|Field/Field_in_season)
# they are account for a sampling design with a time series element and plots

# If we had a designed experiment, we may wish to restrict the permutations so
# that the observations only are permuted within levels of Moisture. Restricted
# permutation is based on the powerful permute package. Function how() can
#be used to define permutation schemes. In the following, we set the levels with
#plots argument:
#  > how <- how(nperm=499, plots = Plots(strata=meta.info$Loc))
# > anova(ord, by="term", permutations = how)

# Do a LDA

# Do a MEM