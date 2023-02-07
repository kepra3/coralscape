# Title: Assessing differentiation of microhabitat by taxa using linear models
# Author: Katharine Prata
# Date: 25/5/22

# Packages
library(tidyverse)
library(lme4)
library(emmeans)
library(ggplot2)
library(vegan)
library(MASS)
library(FSA)

# Functions
combine_metadata <- function(struc.complex, taxa, coordinates) {
  # match
  struc.complex.subset <- subset(struc.complex, sample_name %in% taxa$Individual)
  taxa.subset <- subset(taxa, Individual %in% struc.complex.subset$sample_name)
  coordinates.subset <- subset(coordinates, X %in% struc.complex.subset$sample_name)
  
  # sort
  sorted.taxa.subset <- taxa.subset[order(taxa.subset[,1]),]
  sorted.coordinates.subset <- coordinates.subset[order(coordinates.subset[,1]),]
  sorted.struc.complex.subset <- struc.complex.subset[order(struc.complex.subset[,2]),]
  # check if identiical
  cat("\nChecking taxa and struc complex: ", identical(sorted.taxa.subset[,1], sorted.struc.complex.subset[,2]))
  cat("\nChecking taxa and coordinates: ", identical(sorted.taxa.subset[,1], sorted.coordinates.subset[,1]))
  # Put metadata together
  metadata <- cbind(sorted.taxa.subset, sorted.struc.complex.subset[,-c(1,2)], sorted.coordinates.subset[,-1])
  return(metadata)
}
fix_angles <- function(angles){
  angles_new <- angles
  for (i in 1:length(angles)) {
    if (angles[i] > 90) {
      angles_new[i] <- 180 - angles[i]
    } else {
      angles_new[i] <- angles[i]
    }
  }
  return(angles_new)}
k_boxplot <- function(df, categorical, continuous, lab, colours, facet="No") {
  categorical <- droplevels(categorical)
  present_taxa <- levels(categorical)
  sample.sizes <- NA
  
  for (i in 1:length(present_taxa)) {
    sample.sizes[i] <- length(categorical[categorical == present_taxa[i]])
  }
  
  
  if (facet == "All") {
    p <- ggplot(df, aes(x = categorical, y = continuous, color = categorical)) +
      geom_boxplot(notch = TRUE) +
      geom_point(size = 3, shape = 1) +
      #scale_x_discrete(labels = label) +
      scale_color_manual(name = "Taxa", values = colours,
                         label = present_taxa) +
      theme_classic(base_size = 12) +
      ylab(lab) +
      xlab('Taxa') + 
      facet_wrap(~ Depth + Loc)
  } else if (facet == "Depth") {
    p <- ggplot(df, aes(x = categorical, y = continuous, color = categorical)) +
      geom_boxplot(notch = TRUE) +
      geom_point(size = 3, shape = 1) +
      #scale_x_discrete(labels = label) +
      scale_color_manual(name = "Taxa", values = colours,
                         label = present_taxa) +
      theme_classic(base_size = 12) +
      ylab(lab) +
      xlab('Taxa') + 
      facet_wrap(~ Depth)
  } else if (facet == "Loc") {
    p <- ggplot(df, aes(x = categorical, y = continuous, color = categorical)) +
      geom_boxplot(notch = TRUE) +
      geom_point(size = 3, shape = 1) +
      #scale_x_discrete(labels = label) +
      scale_color_manual(name = "Taxa", values = colours,
                         label = present_taxa) +
      theme_classic(base_size = 12) +
      ylab(lab) +
      xlab('Taxa') + 
      facet_wrap(~ Loc)
  } else {
    label = paste0(present_taxa, "\n n = ", sample.sizes)
    p <- ggplot(df, aes(x = categorical, y = continuous, colour = categorical)) +
      geom_boxplot(notch = TRUE) +
      geom_point(size = 3, shape = 1) +
      scale_x_discrete(labels = label) +
      scale_color_manual(name = "Taxa", values = colours,
                         label = present_taxa) +
      theme_classic(base_size = 12) +
      ylab(lab) +
      xlab('Taxa') }
  return(p)
}

# Import datasets
setwd("~/git/coralscape/results/")
struc.complex <- read.csv("struc_complex_results_envonly.txt", sep = "\t")
taxa <- read.csv("taxa_metadata.csv")
coordinates <- read.csv("annotations.csv")

# Combine
metadata <- combine_metadata(struc.complex, taxa, coordinates)
rm(taxa, coordinates, struc.complex)
write.csv(metadata, "all_metadata_25-5-22.csv", quote = FALSE, row.names = FALSE)

# Remove admixed individuals
metadata <- metadata[!metadata$Taxa == "admixed",]

# Structure ####
str(metadata)
              #AA1        #AA2        #AH1       #AH2      #AH3      #AL1      #AL2
colours <- c("#274e13","#8fce00", "#b06100", "#ffa500", "#ff4e00", "#6f0641", "#7223f0")
metadata[, c(2:5)] <- data.frame(lapply(metadata[, c(2:5)], as.factor))
str(metadata)

# Transformations ####
metadata$environment_rugosity_sqrt <- sqrt(metadata$environment_rugosity)
metadata$ground_elevation_corr <- fix_angles(metadata$ground_elevation)
# Fix values were overhang was greater than colony
metadata$overhang_prop_2 <- NA
for (i in 1:length(metadata$overhang_prop)) {
  if (metadata$overhang_prop[i] > 1) {
    metadata$overhang_prop_2[i] <- 1
  } else {
    metadata$overhang_prop_2[i] <- metadata$overhang_prop[i]
  }
}

# Misassignments ####
to_remove <- NA
# Individuals that assigned to LM but are AC or HU
to_remove <- c("KP0529_LM_WP20", "KP0639_LM_SB20", "KP0399_LM_WP20", "KP0583_LM_WP20")
# Individuals that assigned to AC but are LM (or HU) 
to_remove <- append(to_remove, c("KP0765_AC_WP10", "KP0851_AC_WP10", "KP0108_AC_SB10", "KP0096_AC_SB10"))
# Individual that appears to be AA2, only next one along (771) but assigned to AA1, 590 AA2 but seems like AA1 need to check
to_remove <- append(to_remove, c("KP0769_AC_WP10", "KP0771_AC_WP10", "KP0590_AC_WP20"))
# Funky colony pics, the 650 lamarcki just some adjustment but mostly correct
to_remove <- append(to_remove, c( "KP0650_LM_SB20", "KP0774_AC_WP10", "KP0400_LM_WP20", "KP0350_LM_WP20"))
# Strangely clones with LM
to_remove <- append(to_remove, c("KP0604_AC_SB20"))

metadata <- metadata[!metadata$Individual %in% to_remove,]


# Plots ####
k_boxplot(metadata, metadata$Taxa, metadata$environment_rugosity_sqrt, "Square-root of Environment Rugosity", colours)
k_boxplot(metadata, metadata$Taxa, metadata$environment_rugosity_sqrt, "Square-root of Environment Rugosity", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$environment_rugosity_sqrt, "Square-root of Envirosnment Rugosity", colours, facet = "Depth")

k_boxplot(metadata, metadata$Taxa, metadata$ground_elevation_corr, "Ground Elevation corr", colours)
k_boxplot(metadata, metadata$Taxa, metadata$ground_elevation_corr, "Ground Elevation corr", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$ground_elevation_corr, "Ground Elevation corr", colours, facet = "Depth")

k_boxplot(metadata, metadata$Taxa, metadata$overhang_prop_2, "Overhang proportion", colours)
k_boxplot(metadata, metadata$Taxa, metadata$overhang_prop_2, "Overhang proportion", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$overhang_prop_2, "Overhang proportion", colours, facet = "Depth")

k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours)
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "Depth")


# Loop for saving plots - colours not correct yet for plots that are missing spp.
param_list = c('environment_rugosity_sqrt','outcrop_prop', 'overhang_prop_2', 'ground_elevation_corr')
site_list = c('SB05', 'SB10', 'SB20', 'WP05', 'WP10', 'WP20', "CA05")
depth_list = c('20', '10', '5')

for  (i in site_list) {
  metadata$z_scale[metadata$Site == i] <- scale(metadata$z[metadata$Site == i])
  
}

k_boxplot(metadata, metadata$Taxa, metadata$z_scale, "Local height", colours)
k_boxplot(metadata, metadata$Taxa, metadata$z_scale, "Local height", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$z_scale, "Local height", colours, facet = "Depth")


# Not necessary to run loop
#for (i in param_list) {
#  print("Starting new parameter loop")
#  p <- k_boxplot(metadata, metadata$Taxa, metadata[,i], i, colours)
#  print(p)
#  ggsave(paste0("plots/", i, ".pdf"), p, height = 15, width = 20, units = "cm")
#  p1 <- k_boxplot(metadata, metadata$Taxa, metadata[,i], i, colours, facet = "All")
#  print(p1)
#  ggsave(paste0("plots/", i, "_facet.pdf"), p1, height = 30, width = 20, units = "cm")
#  print("Starting site loop")
#  pd <- k_boxplot(metadata, metadata$Taxa, metadata[,i], i, colours, facet = "Depth")
#  print(pd)
#  ggsave(paste0("plots/", i, "_depth.pdf"), pd, height = 15, width = 20, units = "cm")
#  pl <- k_boxplot(metadata, metadata$Taxa, metadata[,i], i, colours, facet = "Loc")
#  print(pl)
#  ggsave(paste0("plots/", i, "_loc.pdf"), pl, height = 15, width = 20, units = "cm")
#  for (j in site_list) {
#    site <- metadata[metadata$Site == j,]
#    p2 <- k_boxplot(site, site$Taxa, site[,i], paste(i, j), colours)
#    print(p2)
#    ggsave(paste0("plots/", i, "_", j, ".pdf"), p2, height = 15, width = 20, units = "cm")}
#  print("Starting depth loop")
#  for (k in depth_list) {
#    depth <- metadata[metadata$Depth == k,]
#    p3 <- k_boxplot(depth, depth$Taxa, depth[,i], paste(i, k), colours)
#    print(p3)
#    ggsave(paste0("plots/", i, "_", k, ".pdf"), p3, height = 15, width = 20, units = "cm")}
#}

# Proportion of Taxa by Depth ####
table(metadata$Taxa, metadata$Depth)
ggplot(metadata, aes(Depth, fill = Taxa)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = colours) +
  scale_y_continuous(label = scales::percent) +
  theme(axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18, face = 'bold')) +
  ylab("Proportion") +
  theme_classic()

ggplot(metadata, aes(Taxa, fill = Depth)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = blues9[c(3, 6, 9)]) +
  scale_y_continuous(label = scales::percent) +
  theme(axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18, face = 'bold')) +
  ylab("Proportion") +
  theme_classic()

# *** Linear Models ***  ####
## Ground elevation ####
# Normal linear model - not random effects
lm.ground_elevation <- lm(ground_elevation_corr ~ Taxa, data = metadata)
summary(lm.ground_elevation)
TukeyHSD(aov(lm.ground_elevation)) # No significance
#             diff       lwr       upr     p adj
# AL2-AA2 -11.6621075 -24.37905  1.054838 0.0964189
# Is there an effect of Depth
lm.ground_elevation <- lm(ground_elevation_corr ~ Depth, data = metadata)
summary(lm.ground_elevation)
TukeyHSD(aov(lm.ground_elevation))
# No effect of depth

# Random effect of Location
lmer.ground_elevation <- lmer(ground_elevation_corr ~ Taxa + (1|Loc), data = metadata)
summary(lmer.ground_elevation)
emmeans(lmer.ground_elevation, list(pairwise ~ Taxa), adjust = "tukey")
# Is there an effect of Depth - with random effect of location
lmer.ground_elevation.depth <- lmer(ground_elevation_corr ~ Depth + (1|Loc), data = metadata)
summary(lmer.ground_elevation.depth)
emmeans(lmer.ground_elevation.depth, list(pairwise ~ Depth), adjust = "tukey")
# No effect of depth

# Interaction b/w Taxa and Depth
lmer2.ground_elevation <- lmer(ground_elevation_corr ~ Taxa*Depth + (1|Loc), data = metadata)
summary(lmer2.ground_elevation)
emmeans(lmer2.ground_elevation, list(pairwise ~ Taxa:Depth), adjust = "tukey") # No interactions

lmer3.ground_elevation <- lmer(ground_elevation_corr ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer3.ground_elevation)
emmeans(lmer3.ground_elevation, list(pairwise ~ Taxa), adjust = "tukey")

lmer4.ground_elevation <- lmer(ground_elevation_corr ~ Taxa + (Taxa|Loc/Depth), data = metadata)
summary(lmer4.ground_elevation)
emmeans(lmer4.ground_elevation, list(pairwise ~ Taxa), adjust = "tukey") # No interactions
# NaNs

# tests with difference variance...
kruskal.test(metadata$ground_elevation_corr, metadata$Taxa)
boxplot(metadata$ground_elevation_corr ~ metadata$Taxa)
# Not significant
#dunnTest(metadata$ground_elevation ~ metadata$Taxa, 
#         method = 'bonferroni')

## Environment rugosity ####
lm.environment_rugosity <- lm(environment_rugosity_sqrt ~ Taxa, data = metadata)
summary(lm.environment_rugosity)
TukeyHSD(aov(lm.environment_rugosity))
# AL2-AA2 -0.123497839 -0.21245816 -0.03453752 0.0009383
# AH2-AA2 -0.158586302 -0.30013033 -0.01704227 0.0169592
# AH3-AA2 -0.185234211 -0.35395456 -0.01651386 0.0209918
boxplot(metadata$environment_rugosity_sqrt ~ metadata$Taxa)
# Most different to AA2 (with high rugosity)
plot(lm.environment_rugosity)
# Depth
lm.environment_rugosity <- lm(environment_rugosity_sqrt ~ Depth, data = metadata)
summary(lm.environment_rugosity)
TukeyHSD(aov(lm.environment_rugosity))
#               diff         lwr        upr     p adj
#10-5  0.1392920042  0.06610708 0.21247693 0.0000301
#20-5  0.1395009794  0.07225033 0.20675163 0.0000048
#20-10 0.0002089752 -0.04811260 0.04853055 0.9999429
# Low rugosity at 5m

# Location
lm.environment_rugosity <- lm(environment_rugosity_sqrt ~ Loc, data = metadata)
summary(lm.environment_rugosity)
TukeyHSD(aov(lm.environment_rugosity))
#diff         lwr        upr     p adj
#SB-CA 0.10231454 -0.01349110 0.21812018 0.0955674
#WP-CA 0.12135850  0.00900138 0.23371562 0.0306714
#WP-SB 0.01904396 -0.02830004 0.06638797 0.6110456
boxplot(metadata$environment_rugosity_sqrt ~ metadata$Loc)
# Low rugosity at CA (but indivduals at CA are only from 5m ...)

# Including Loc/Depth as random effect
lmer.environment_rugosity <- lmer(environment_rugosity_sqrt ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer.environment_rugosity)
emmeans(lmer.environment_rugosity, list(pairwise ~ Taxa), adjust = "tukey")
#  AA2 - AL2   0.1352 0.0305 345   4.431  0.0003

lmer.environment_rugosity.depth <- lmer(environment_rugosity_sqrt ~ Depth + (1|Loc), data = metadata)
summary(lmer.environment_rugosity.depth)
emmeans(lmer.environment_rugosity.depth, list(pairwise ~ Depth), adjust = "tukey")

# Interaction
lmer2.environment_rugosity <- lmer(environment_rugosity_sqrt ~ Taxa * Depth + (1|Loc), data = metadata)
summary(lmer2.environment_rugosity)
emmeans(lmer2.environment_rugosity, list(pairwise ~ Taxa:Depth), adjust = "tukey")
# Doesn't neccessary have to do with the Taxa just that environment rugosity changes with depth
env_rug_depth.taxa.int <- emmeans(lmer2.environment_rugosity, list(pairwise ~ Taxa:Depth), adjust = "tukey")
options(max.print = 10000)
env_rug_depth.taxa.int$`pairwise differences of Taxa, Depth`
# AA2 Depth20 - AL2 Depth20  0.12710 0.0329 328   3.859  0.0210
# AA2 Depth5 - AA2 Depth10  -0.08956 0.0480 274  -1.868  0.9536

# Conclusions: Rugosity changes with depth, 5 m less rugose and 10 and 20 m equally rugose.
# At 20 m, AA2 occupies a more rugose environment than AL2
# But the colony itself is included in this measure so probably just a different in colony morphology

## Outcrop proportion ####
lm.outcrop_prop <- lm(outcrop_prop ~ Taxa, data = metadata)
summary(lm.outcrop_prop)
TukeyHSD(aov(lm.outcrop_prop))
# AL1-AA2 -0.094684990 -0.18572171 -0.003648269 0.0354633
# AL2-AA2 -0.136333539 -0.19754460 -0.075122478 0.0000000
# AL2-AA1 -0.113673334 -0.20614510 -0.021201566 0.0056501
boxplot(metadata$outcrop_prop ~ metadata$Taxa * metadata$Depth)

lm.outcrop.depth <- lm(outcrop_prop ~ Depth, data = metadata)
summary(lm.outcrop.depth)
TukeyHSD(aov(lm.outcrop.depth))
# 20-10  0.03723189  0.002388313 0.07207547 0.0329851
boxplot(metadata$outcrop_prop ~ metadata$Depth)

lmer.outcrop_prop <- lmer(outcrop_prop ~ Taxa*Depth + (1|Loc), data = metadata)
summary(lmer.outcrop_prop)
emmeans(lmer.outcrop_prop, list(pairwise ~ Taxa:Depth), adjust = "tukey")
# AA2 Depth10 - AA2 Depth20 -0.059731 0.0158 337  -3.783  0.0271
# AA1 Depth20 - AL2 Depth20  0.108189 0.0308 337   3.511  0.0656
# AH1 Depth10 - AA2 Depth20 -0.141223 0.0347 337  -4.066  0.0098
# AA2 Depth10 - AL2 Depth20  0.096426 0.0229 337   4.205  0.0057

lmer.outcrop_prop.depth <- lmer(outcrop_prop ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer.outcrop_prop.depth)
emmeans(lmer.outcrop_prop.depth, list(pairwise ~Taxa), adjust = "tukey")
# AH2 - AL2  0.13015 0.0429 135   3.034  0.0447
# AA2 - AL1  0.10976 0.0305 345   3.602  0.0066
# AA2 - AL2  0.14756 0.0207 345   7.118  <.0001
# AA1 - AL2  0.09952 0.0307 343   3.244  0.0218


glmer.outcrop.interaction <- glmer(outcrop_prop ~ Taxa + (1|Loc/Depth), family = 'binomial', data = metadata)
emmeans(glmer.outcrop.interaction, list(pairwise~Taxa))
# AA2 - AH1   1.8683 0.435 Inf   4.290  0.0004
# AA2 - AL2   1.8499 0.379 Inf   4.875  <.0001
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "Depth")

# Conclusions: probably most interesting one, AA2 seem to change between 10-20,
#& AH1 it's relation to AA2
# AA1 slightly lower (no different to AL1 but diff to AL2), AL1 and AL2 lowest.
# but it is also to do with colony morphology (i.e., being top heavy?)
plot(metadata$outcrop_prop ~ metadata$z)

## Overhang proportion ####
lm.overhang <- lm(overhang_prop_2 ~ Taxa, data = metadata)
summary(lm.overhang)
anova(lm.overhang)
TukeyHSD(aov(lm.overhang))
# AH1-AA2  0.13197544 -0.003130884 0.26708176 0.0605199
lm.overhang.depth <- lm(overhang_prop_2 ~ Depth, data = metadata)
summary(lm.overhang.depth)
anova(lm.overhang.depth)
TukeyHSD(aov(lm.overhang.depth))

lm.overhang.loc <- lm(overhang_prop_2 ~ Loc, data = metadata)
summary(lm.overhang.loc)
anova(lm.overhang.loc)
TukeyHSD(aov(lm.overhang.loc))

lmer.overhang <- lmer(overhang_prop_2 ~ Taxa + (1|Loc), data = metadata)
summary(lmer.overhang)
emmeans(lmer.overhang, list(pairwise~Taxa), adjust = "tukey")

lmer2.overhang <- lmer(overhang_prop_2 ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer2.overhang)
emmeans(lmer2.overhang, list(pairwise~Taxa), adjust = "tukey")
# environment radius probably too small to detect most overhangs!!!
k_boxplot(metadata, metadata$Taxa, metadata$overhang_prop_2, "Overhang proportion", colours, facet = "Depth")

glmer.overhang <- glmer(overhang_prop_2 ~ Taxa + (1|Loc/Depth), family = 'binomial', data = metadata)
emmeans(glmer.overhang, list(pairwise~Taxa))

# Conclusions: AH1 likes to be cryptic at 10m. AL1 and AL2 probably have higher overhangs when increasing radius.
# General patterns make sense

## Z - local depth... ####
# Linear models
lm.z <- lm(z_scale ~ Taxa, data = metadata)
summary(lm.z)
anova(lm.z)
TukeyHSD(aov(lm.z))
# AL2-AA2 -0.88278968 -1.3893180 -0.37626134 0.0000083
# AL2-AH1 -0.73032152 -1.4580838 -0.00255929 0.0485478
# AL2-AH2 -1.10333055 -2.0183002 -0.18836086 0.0072233
# AL2-AH3 -1.47983769 -2.5336550 -0.42602034 0.0007736

lm.z.depth <- lm(z_scale ~ Depth, data = metadata)
summary(lm.z.depth)
anova(lm.z.depth)
TukeyHSD(aov(lm.z.depth))

lm.z.loc <- lm(z_scale ~ Loc, data = metadata)
summary(lm.z.loc)
anova(lm.z.loc)
TukeyHSD(aov(lm.z.loc))
# Shouldn't be any differences between Depth and Location (alone) because they are all scaled

lmer.z <- lmer(z_scale ~ Taxa + (1|Loc), data = metadata)
summary(lmer.z)
emmeans(lmer.z, list(pairwise~Taxa), adjust = "tukey")

lmer2.z <- lmer(z_scale ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer2.z)
emmeans(lmer2.z, list(pairwise~Taxa), adjust = "tukey")

k_boxplot(metadata, metadata$Taxa, metadata$z_scale, "Local height", colours, facet = "Depth")
# a bit different to outcrop proportions
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "Depth")
