# Title: Niche differences
# Author: Katharine Prata
# Date created: 17/01/21

# Packages
library(vcfR)
library(adegenet)
library(ggplot2)
library(lme4)
library(emmeans)

# Functions
combine_metadata <- function(struc.complex, clusters, coordinates) {
  # match
  struc.complex.subset <- subset(struc.complex, X %in% clusters$Individual)
  clusters.subset <- subset(clusters, Individual %in% struc.complex.subset$X)
  coordinates.subset <- subset(coordinates, X %in% struc.complex.subset$X)
  
  # sort
  sorted.clusters.subset <- clusters.subset[order(clusters.subset[,1]),]
  sorted.coordinates.subset <- coordinates.subset[order(coordinates.subset[,1]),]
  sorted.struc.complex.subset <- struc.complex.subset[order(struc.complex.subset[,1]),]
  # check if identiical
  cat("\nChecking clusters and struc complex: ", identical(sorted.clusters.subset[,1], sorted.struc.complex.subset[,1]))
  cat("\nChecking clusters and coordinates: ", identical(sorted.clusters.subset[,1], sorted.coordinates.subset[,1]))
  # Put metadata together
  metadata <- sorted.clusters.subset
  metadata[,20:25] <- sorted.struc.complex.subset[,2:7]
  metadata[, 26:28] <- sorted.coordinates.subset[,2:4]
  return(metadata)
}

overhang_barplot <- function(metadata, factor) {
  overhang <- table(metadata$overhang, factor) %>% matrix(nrow = 2, ncol = length(levels(factor)))
  colnames(overhang) <- levels(factor)
  rownames(overhang) <- levels(metadata$overhang)
  barplot(overhang, main = "Presence of overhang", xlab = "Taxa",
          col = c("Red", "Green"))
  legend("topright", c("No", "Yes"), fill = c("Red", "Green"))
}

k_boxplot <- function(df, cat, cont) {
  
  levels <- levels(cat)
  sample.sizes <- NA
  
  for (i in 1:length(levels)) {
    sample.sizes[i] <- length(cat[cat == levels[i]])
  }
  
  label = paste0(levels, "\n n = ", sample.sizes)
  
  p <- ggplot(df, aes(x = cat, y = cont)) +
    geom_boxplot(notch = TRUE) +
    geom_point(size = 3, shape = 1) +
    scale_x_discrete(labels = label) +
    theme_classic(base_size = 12)
  return(p)
}

# Import datasets
struc.complex <- read.csv("~/git/coralscape/results/sample_metadata.csv")
clusters <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/all-aga_1d_nc_20_6.csv")
coordinates <- read.csv("~/git/coralscape/results/annotations.csv")

# Organise data
metadata <- combine_metadata(struc.complex, clusters, coordinates)

# Clean working space
rm(struc.complex)
rm(clusters)
rm(coordinates)

# Data structure
str(metadata)
metadata[, c(2:6, 13)] <- data.frame(lapply(metadata[, c(2:6, 13)], as.factor))
levels(metadata$Clusters) <- c("AA1", "AL2", "AH1", "AA2", "AL1","AH2")
str(metadata)

# choosing steepest angle
metadata$theta <- NA
for (row in 1:length(metadata$xz)) {
  if (is.na(metadata$xz[row])) {
    next
  } else if (is.na(metadata$yz[row])) {
    next
  } else if (metadata$xz[row] > metadata$yz[row]) {
    metadata$theta[row] <- metadata$xz[row]
  } else {
    metadata$theta[row] <- metadata$yz[row]
  }
}

# Environmental data ####
# Overhang
boxplot(metadata$overhang_prop ~ metadata$Clusters)

# Theta
k_boxplot(metadata, metadata$Species, metadata$theta) + ggtitle('Theta')
k_boxplot(metadata, metadata$Clusters, metadata$theta) + ggtitle('Theta')
# maybe colony rugosity would be good to include here!

# Relative depth proportion species
k_boxplot(metadata, metadata$Species, metadata$outcrop_prop) + ggtitle('Outcrop position')
k_boxplot(metadata, metadata$Clusters, metadata$outcrop_prop) + ggtitle('Outcrop position')

# Z
k_boxplot(metadata, metadata$Species, metadata$z) + ggtitle('Relative depth') # doesn't mean anything for now
k_boxplot(metadata, metadata$Clusters, metadata$z) + ggtitle('Relative depth')

#### Linear Models ####
# 1. overhang #### will become overhang %
metadata$overhang.bin <- ifelse(metadata$overhang == "Yes", 1, 0)
glmer1.overhang <- glmer(overhang.bin ~ Clusters + (1|Loc/Depth), data = metadata, family = binomial) # won't work
# boundary (singular) fit: see ?isSingular - means the random effects don't use too much variation
summary(glmer1.overhang)
#Warning messages:
#  1: In vcov.merMod(object, use.hessian = use.hessian) :
#  variance-covariance matrix computed from finite-difference Hessian is
#not positive definite or contains NA values: falling back to var-cov estimated from RX
#2: In vcov.merMod(object, correlation = correlation, sigm = sig) :
#  variance-covariance matrix computed from finite-difference Hessian is
#not positive definite or contains NA values: falling back to var-cov estimated from RX
# ClustersAL1  2.807e+00  6.065e-01   4.628 3.69e-06 ***
# ClustersAL2  2.750e+00  4.819e-01   5.706 1.16e-08 ***

emmeans(glmer1.overhang, list(pairwise ~ Clusters), adjust = "tukey")
#AA1 - AL2  -2.7496        0 Inf -5.706  <.0001 ***
#AL2 - AH1   3.3101        1 Inf  3.088  0.0246 *
#AH1 - AL1  -3.3673        1 Inf -2.971  0.0352 *

# probability space with confidence interval?

# 2. theta ####

shapiro.test(metadata$theta) # not normal

lm.theta <- lm(theta ~ Clusters, data = metadata)
summary(lm.theta) # all diff to AA1 except AH2 :D
summary(aov(lm.theta))
TukeyHSD(aov(lm.theta))
#AL2-AA1  5.8931198   0.01654314 11.769696 0.0489002 *
#AH1-AA1  9.7738198   3.50915167 16.038488 0.0001624 ***
#AA2-AA1  8.3621034   0.10943483 16.614772 0.0449977 *
#AL1-AA1  9.6290353   1.37636675 17.881704 0.0118319 *



# mixed model
lmer1.theta <- lmer(theta ~ Clusters + (1|Loc/Depth), data = metadata)
# boundary (singular) fit: see ?isSingular - means that the random effects don't use too much variation
summary(lmer1.theta)
# extract coefficients
coefs <- data.frame(coef(summary(lmer1.theta)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
# post-hoc
emmeans(lmer1.theta, list(pairwise ~ Clusters), adjust = "tukey")
#AA1 - AL2   -5.893 2.15 237 -2.737  0.0718 `
#AA1 - AH1   -9.774 2.41 164 -4.054  0.0011 **
#AA1 - AA2   -8.362 3.02 232 -2.771  0.0659 `
#AA1 - AL1   -9.629 2.94 249 -3.274  0.0152 * 

# mixed model 2
lmer2.theta <- lmer(theta ~ Clusters + (Clusters|Loc/Depth), data = metadata)
# boundary (singular) fit: see ?isSingular - means that the random effects don't use too much variation
summary(lmer2.theta)
# extract coefficients
coefs <- data.frame(coef(summary(lmer2.theta)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
# post-hoc
emmeans(lmer2.theta, list(pairwise ~ Clusters), adjust = "tukey")
# Warning message:
#  In ptukey(sqrt(2) * abst, fam.size, zapsmall(df), lower.tail = FALSE) :
#  NaNs produced

# including colony size ####
lmer3.theta <- lmer(theta ~ Clusters * range + (1|Loc/Depth), data = metadata)
summary(lmer3.theta)
coefs <- data.frame(coef(summary(lmer3.theta)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
emmeans(lmer3.theta, list(pairwise ~ Clusters), adjust = "tukey")


# including overhang ####
lmer4.theta <- lmer(theta ~ Clusters * overhang + (1|Loc/Depth), data = metadata)
summary(lmer4.theta)
coefs <- data.frame(coef(summary(lmer4.theta)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
emmeans(lmer4.theta, list(pairwise ~ Clusters), adjust = "tukey")
# maybe issue because some species have no overhang presence

# 3. outcrop proportion ####
hist(metadata[metadata$Clusters == "AA1",]$prop)
hist(metadata[metadata$Clusters == "AL2",]$prop)
hist(metadata[metadata$Clusters == "AH1",]$prop)
hist(metadata[metadata$Clusters == "AA2",]$prop)
hist(metadata[metadata$Clusters == "AL1",]$prop)
hist(metadata[metadata$Clusters == "AH2",]$prop)
# sort of normally distributed
boxplot(prop ~ Clusters, data = metadata)

lm.prop <- lm(prop ~ Clusters, data = metadata)
summary(lm.prop) # only AL1 and AL2 diff to AA1
summary(aov(lm.prop))
TukeyHSD(aov(lm.prop))
# AL2-AA1 -0.122970975 -0.180950729 -0.064991220 0.0000001 ***
# AL1-AA1 -0.091952545 -0.173375409 -0.010529682 0.0167064 *
# AH1-AL2  0.118848284  0.041629615  0.196066952 0.0002101 ***
# AA2-AL2  0.098519324  0.004859928  0.192178721 0.0327950 *


lmer1.prop <- lmer(prop ~ Clusters + (1|Loc/Depth), data = metadata )
summary(lmer1.prop)
# extract coefficients
coefs <- data.frame(coef(summary(lmer1.prop)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
# post-hoc
emmeans(lmer1.prop, list(pairwise ~ Clusters), adjust = "tukey")
# AA1 - AL2  0.122882 0.0209 255  5.879  <.0001 ***
# AA1 - AL1  0.092556 0.0288 256  3.213  0.0184 *
# AL2 - AH1 -0.122257 0.0307 173 -3.983  0.0014 **

# 4. rugosity