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
#struc.complex <- read.csv("env_0.2/all_struc_complex.txt", sep = "\t")
#struc.complex <- struc.complex[struc.complex$plot_name != "plot_name",]
#taxa <- read.csv("taxa_metadata.csv")
#coordinates <- read.csv("all_annotations_X_DEPTH_parallel.txt", sep = "\t")

# Combine
#metadata <- combine_metadata(struc.complex, taxa, coordinates)
#rm(taxa, coordinates, struc.complex)
#write.csv(metadata, "all_metadata_14-6-22.csv", quote = FALSE, row.names = FALSE)
metadata <- read.csv("all_metadata_14-6-22.csv")
# Remove admixed individuals
metadata <- metadata[!metadata$Taxa == "admixed",]

# Structure ####
str(metadata)
#AA1        #AA2        #AH1       #AH2      #AH3      #AL1      #AL2
colours <- c("#274e13","#8fce00", "#b06100", "#ffa500", "#ff4e00", "#6f0641", "#7223f0")
metadata$Depth[metadata$Depth == "10"] <- "10-12"
metadata$Depth[metadata$Depth == "12"] <- "10-12"
metadata$Depth <- factor(metadata$Depth, levels = c("5", "10-12", "20"))
metadata[, c(2:6)] <- data.frame(lapply(metadata[, c(2:6)], as.factor))
metadata$colony_points <- as.integer(metadata$colony_points)
metadata[, 8:15] <- data.frame(lapply(metadata[, 8:15], as.numeric))
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
site_list = c('SB05', 'SB10', 'SB20', 'WP05', 'WP10', 'WP20', "CA05", "CA10", "CA20", "SQ12", "SQ20")
depth_list = c('20', '10-12', '5')

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
chisq.test(table(metadata$Taxa, metadata$Depth))

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

# sizes ####

mean(metadata$colony_range[metadata$Taxa == "AA1"])
sd(metadata$colony_range[metadata$Taxa == "AA1"])
mean(metadata$colony_range[metadata$Taxa == "AA2"])
sd(metadata$colony_range[metadata$Taxa == "AA2"])

mean(metadata$colony_range[metadata$Taxa == "AH1"])
sd(metadata$colony_range[metadata$Taxa == "AH1"])
mean(metadata$colony_range[metadata$Taxa == "AH2"])
sd(metadata$colony_range[metadata$Taxa == "AH2"])
mean(metadata$colony_range[metadata$Taxa == "AH3"])
sd(metadata$colony_range[metadata$Taxa == "AH3"])

mean(metadata$colony_range[metadata$Taxa == "AL1"])
sd(metadata$colony_range[metadata$Taxa == "AL1"])
mean(metadata$colony_range[metadata$Taxa == "AL2"])
sd(metadata$colony_range[metadata$Taxa == "AL2"])


# colony rugosity ####
mean(metadata$colony_rugosity[metadata$Taxa == "AA1"])
sd(metadata$colony_rugosity[metadata$Taxa == "AA1"])
mean(metadata$colony_rugosity[metadata$Taxa == "AA2"])
sd(metadata$colony_rugosity[metadata$Taxa == "AA2"])

mean(metadata$colony_rugosity[metadata$Taxa == "AH1"])
sd(metadata$colony_rugosity[metadata$Taxa == "AH1"])
mean(metadata$colony_rugosity[metadata$Taxa == "AH2"])
sd(metadata$colony_rugosity[metadata$Taxa == "AH2"])
mean(metadata$colony_rugosity[metadata$Taxa == "AH3"])
sd(metadata$colony_rugosity[metadata$Taxa == "AH3"])

mean(metadata$colony_rugosity[metadata$Taxa == "AL1"])
sd(metadata$colony_rugosity[metadata$Taxa == "AL1"])
mean(metadata$colony_rugosity[metadata$Taxa == "AL2"])
sd(metadata$colony_rugosity[metadata$Taxa == "AL2"])

# *** Linear Models ***  ####
## Ground elevation ####
k_boxplot(metadata, metadata$Taxa, metadata$ground_elevation_corr, "Ground Elevation corr", colours)
k_boxplot(metadata, metadata$Taxa, metadata$ground_elevation_corr, "Ground Elevation corr", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$ground_elevation_corr, "Ground Elevation corr", colours, facet = "Depth")
ggsave("ground_elevation_by_site.pdf", height = 11, width = 20, units = "cm")
# Normal linear model - not random effects
lm.ground_elevation <- lm(ground_elevation_corr ~ Taxa, data = metadata)
summary(lm.ground_elevation)
TukeyHSD(aov(lm.ground_elevation))
#             diff       lwr       upr     p adj
# AL2-AH1 -15.097466614 -30.620823  0.4258899 0.0627055

# Is there an effect of Depth
lm.ground_elevation <- lm(ground_elevation_corr ~ Depth, data = metadata)
summary(lm.ground_elevation)
TukeyHSD(aov(lm.ground_elevation))
# Yes - higher slopes at 20 m
plot(metadata$ground_elevation_corr ~ metadata$Depth)

# Random effect of Location/Depth
lmer.ground_elevation <- lmer(ground_elevation_corr ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer.ground_elevation)
emmeans(lmer.ground_elevation, list(pairwise ~ Taxa), adjust = "tukey")
# AH1 - AL2 16.54095 5.47 414   3.022  0.0423
# AA2 - AL2  9.69854 3.44 494   2.816  0.0744
# Is there an effect of Depth - with random effect of location
lmer.ground_elevation.depth <- lmer(ground_elevation_corr ~ Depth + (1|Loc/Depth), data = metadata)
summary(lmer.ground_elevation.depth)
emmeans(lmer.ground_elevation.depth, list(pairwise ~ Depth), adjust = "tukey")
# No effect of depth

# Interaction b/w Taxa and Depth
lmer2.ground_elevation <- lmer(ground_elevation_corr ~ Taxa*Depth + (1|Loc/Depth), data = metadata)
summary(lmer2.ground_elevation)
emmeans(lmer2.ground_elevation, list(pairwise ~ Taxa:Depth), adjust = "tukey") # No interactions

# tests with difference variance...
kruskal.test(metadata$ground_elevation_corr, metadata$Taxa)
boxplot(metadata$ground_elevation_corr ~ metadata$Taxa)
# Significant
dunnTest(metadata$ground_elevation_corr ~ metadata$Taxa, 
         method = 'bonferroni')

# conclusion: AH1 sits at higher angles than AL2 (not dependent on depth...)

## Environment rugosity ####
k_boxplot(metadata, metadata$Taxa, metadata$environment_rugosity_sqrt, "Square-root of Environment Rugosity", colours)
k_boxplot(metadata, metadata$Taxa, metadata$environment_rugosity_sqrt, "Square-root of Environment Rugosity", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$environment_rugosity_sqrt, "Square-root of Envirosnment Rugosity", colours, facet = "Depth")
ggsave("environment_rugosity_by_depth.pdf", height = 11, width = 20, units = "cm")

lm.environment_rugosity <- lm(environment_rugosity_sqrt ~ Taxa, data = metadata)
summary(lm.environment_rugosity)
TukeyHSD(aov(lm.environment_rugosity))
# AH2-AA2 -0.152678217 -0.28941136 -0.01594508 0.0174810
# just because more AH2 at 5m
boxplot(metadata$environment_rugosity_sqrt ~ metadata$Taxa)

# Depth
lm.environment_rugosity <- lm(environment_rugosity_sqrt ~ Depth, data = metadata)
summary(lm.environment_rugosity)
TukeyHSD(aov(lm.environment_rugosity))
#                diff         lwr        upr     p adj
#10-12-5   0.16217150  0.10294501 0.22139800 0.0000000
#20-5      0.15072728  0.09473468 0.20671988 0.0000000
#20-10-12 -0.01144422 -0.04554902 0.02266057 0.7100988
# Low rugosity at 5m
boxplot(metadata$environment_rugosity_sqrt ~ metadata$Depth)

# Location
lm.environment_rugosity <- lm(environment_rugosity_sqrt ~ Loc, data = metadata)
summary(lm.environment_rugosity)
TukeyHSD(aov(lm.environment_rugosity))
#diff         lwr        upr     p adj
#SQ-CA -0.08008961 -0.141265724 -0.018913491 0.0044154
#SQ-SB -0.06337168 -0.122415944 -0.004327416 0.0298618

boxplot(metadata$environment_rugosity_sqrt ~ metadata$Loc)
# Low rugosity at SQ

# Including Loc/Depth as random effect
lmer.environment_rugosity <- lmer(environment_rugosity_sqrt ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer.environment_rugosity)
emmeans(lmer.environment_rugosity, list(pairwise ~ Taxa), adjust = "tukey")
# AA2 - AL1  0.099657 0.0295 492   3.377  0.0139
# AA2 - AL2  0.076752 0.0213 493   3.611  0.0062
boxplot(metadata$environment_rugosity_sqrt ~ metadata$Taxa)
# AL perfer lower rugosity than AA2

lmer.environment_rugosity.depth <- lmer(environment_rugosity_sqrt ~ Depth + (1|Loc), data = metadata)
summary(lmer.environment_rugosity.depth)
emmeans(lmer.environment_rugosity.depth, list(pairwise ~ Depth), adjust = "tukey")

# Interaction
lmer2.environment_rugosity <- lmer(environment_rugosity_sqrt ~ Taxa * Depth + (1|Loc/Depth), data = metadata)
summary(lmer2.environment_rugosity)
options(max.print = 10000)
emmeans(lmer2.environment_rugosity, list(pairwise ~ Taxa:Depth), adjust = "tukey")
env_rug_depth.taxa.int <- emmeans(lmer2.environment_rugosity, list(pairwise ~ Taxa:Depth), adjust = "tukey")
env_rug_depth.taxa.int$`pairwise differences of Taxa, Depth`
# No within depth differences.... 5 and (10-12)/20 differences

# Conclusions: Rugosity changes with depth, 5 m less rugose and 10 and 20 m equally rugose.
# AL prefer lower rugosity than AA2?

## Outcrop proportion ####
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours)
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "Depth")
ggsave("outcrop_proportion_by_depth.pdf", height = 11, width = 20, units = "cm")


lm.outcrop_prop <- lm(outcrop_prop ~ Taxa, data = metadata)
summary(lm.outcrop_prop)
TukeyHSD(aov(lm.outcrop_prop))
# AL2-AH1 -0.108266012 -0.19092268 -0.02560934 0.0022842
# AL2-AA2 -0.141935225 -0.19582844 -0.08804201 0.0000000
# AL1-AA2 -0.093074357 -0.16891100 -0.01723771 0.0056865
# AL2-AA1 -0.128126974 -0.20838319 -0.04787076 0.0000609
# AL2-AL1 -0.048860868 -0.13778528  0.04006355 0.6651116

lm.outcrop.depth <- lm(outcrop_prop ~ Depth, data = metadata)
summary(lm.outcrop.depth)
TukeyHSD(aov(lm.outcrop.depth))
# No differences
boxplot(metadata$outcrop_prop ~ metadata$Depth)

lmer.outcrop_prop.depth <- lmer(outcrop_prop ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer.outcrop_prop.depth)
emmeans(lmer.outcrop_prop.depth, list(pairwise ~Taxa), adjust = "tukey")
# AH1 - AL2  0.10476 0.0294 439   3.561  0.0074
# AA1 - AL2  0.11195 0.0270 496   4.144  0.0008
# AA2 - AL1  0.09831 0.0255 496   3.853  0.0025
# AA2 - AL2  0.14299 0.0184 497   7.787  <.0001


lmer.outcrop_prop <- lmer(outcrop_prop ~ Taxa*Depth + (1|Loc/Depth), data = metadata)
summary(lmer.outcrop_prop)
anova(lmer.outcrop_prop)
emmeans(lmer.outcrop_prop, list(pairwise ~ Taxa:Depth), adjust = "tukey")
# AH1 20 - AL2 20            0.118314 0.0486 484.42   2.434  0.6597 - not within depth
# AA1 20 - AL2 20            0.116171 0.0280 487.12   4.147  0.0068
# AA2 20 - AL1 20            0.101419 0.0262 483.32   3.864  0.0197
# AA2 20 - AL2 20            0.152363 0.0199 485.97   7.667  <.0001



glmer.outcrop.interaction <- glmer(outcrop_prop ~ Taxa + (1|Loc/Depth), family = 'binomial', data = metadata)
emmeans(glmer.outcrop.interaction, list(pairwise~Taxa))
#AL1 - AL2   1.4759 0.528 Inf   2.796  0.0763  intersting...
#AH1 - AL2   1.7944 0.511 Inf   3.509  0.0082
#AA2 - AL2   2.0808 0.329 Inf   6.324  <.0001
#AA1 - AL2   2.0950 0.522 Inf   4.014  0.0012

plot(metadata$outcrop_prop ~ metadata$z_scale) # maybe slight correlation...
summary(lm(metadata$outcrop_prop ~ metadata$z_scale)) # slight correlation!

# AL1/AL2 - bottom of outcrop

## Overhang proportion ####
k_boxplot(metadata, metadata$Taxa, metadata$overhang_prop_2, "Overhang proportion", colours)
k_boxplot(metadata, metadata$Taxa, metadata$overhang_prop_2, "Overhang proportion", colours, facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$overhang_prop_2, "Overhang proportion", colours, facet = "Depth")
ggsave("overhang_proportion_by_depth.pdf", height = 11, width = 20, units = "cm")

lm.overhang <- lm(overhang_prop_2 ~ Taxa, data = metadata)
summary(lm.overhang)
anova(lm.overhang)
TukeyHSD(aov(lm.overhang))
# AL2-AA2  0.1365622032  0.018032558 0.2550918 0.0123497
# AH1-AA2  0.1504283687  0.000038186 0.3008186 0.0498934
lm.overhang.depth <- lm(overhang_prop_2 ~ Depth, data = metadata)
summary(lm.overhang.depth)
anova(lm.overhang.depth)
TukeyHSD(aov(lm.overhang.depth))
# No difference with depth

lm.overhang.loc <- lm(overhang_prop_2 ~ Loc, data = metadata)
summary(lm.overhang.loc)
anova(lm.overhang.loc)
TukeyHSD(aov(lm.overhang.loc))
# SQ-SB -0.11173373 -0.216735480 -0.006731982 0.0319196
plot(metadata$overhang_prop_2 ~ metadata$Loc) # less overhangs at SQ

lmer.overhang <- lmer(overhang_prop_2 ~ Taxa + (1|Loc), data = metadata)
summary(lmer.overhang)
emmeans(lmer.overhang, list(pairwise~Taxa), adjust = "tukey")
# AA2 - AH1 -0.14949 0.0512 497  -2.921  0.0559
# AA2 - AL2 -0.13456 0.0401 499  -3.352  0.0150

lmer2.overhang <- lmer(overhang_prop_2 ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer2.overhang)
emmeans(lmer2.overhang, list(pairwise~Taxa), adjust = "tukey")
# AA2 - AL2 -0.13384 0.0410 466  -3.265  0.0200
# AA2 - AH1 -0.14427 0.0525 406  -2.748  0.0893

lmer3.overhang <- lmer(overhang_prop_2 ~ Taxa * Depth + (1|Loc/Depth), data = metadata)
summary(lmer3.overhang)
emmeans(lmer3.overhang, list(pairwise~Taxa:Depth), adjust = "tukey")
# AA2 20 - AL2 20           -0.15638 0.0444 489.81  -3.520  0.0621


glmer.overhang <- glmer(overhang_prop_2 ~ Taxa + (1|Loc/Depth), family = 'binomial', data = metadata)
emmeans(glmer.overhang, list(pairwise~Taxa))
# AA2 - AL2  -1.0348 0.376 Inf  -2.750  0.0863

# Conclusions: AL2 more overhangs, maybe AH1 as well than AA2.

## Z - local depth... ####
# Linear models
k_boxplot(metadata, metadata$Taxa, metadata$z_scale, "Local height", colours, facet = "Depth")
ggsave("local_height_by_depth.pdf", height = 11, width = 20, units = "cm")

lm.z <- lm(z_scale ~ Taxa, data = metadata)
summary(lm.z)
anova(lm.z)
TukeyHSD(aov(lm.z))
# AL2-AA2 -0.798619423 -1.2248412 -0.37239764 0.0000010
# AL2-AH2 -0.975539554 -1.9200290 -0.03105012 0.0377154
# AL2-AH3 -1.276886210 -2.1884294 -0.36534302 0.0007790

lm.z.depth <- lm(z_scale ~ Depth, data = metadata)
summary(lm.z.depth)
anova(lm.z.depth)
TukeyHSD(aov(lm.z.depth))


lm.z.loc <- lm(z_scale ~ Loc, data = metadata)
summary(lm.z.loc)
anova(lm.z.loc)
TukeyHSD(aov(lm.z.loc))

lmer.z <- lmer(z_scale ~ Taxa + (1|Loc), data = metadata)
summary(lmer.z)
emmeans(lmer.z, list(pairwise~Taxa), adjust = "tukey")
# AH3 - AL2  1.27689 0.312 494   4.095  0.0010
# AH2 - AL2  0.97554 0.320 499   3.045  0.0392
# AA2 - AL2  0.79862 0.145 499   5.524  <.0001

lmer2.z <- lmer(z_scale ~ Taxa + (1|Loc/Depth), data = metadata)
summary(lmer2.z)
emmeans(lmer2.z, list(pairwise~Taxa), adjust = "tukey")
# AA1 - AA2 -0.32073 0.182 378  -1.759  0.5764
# AA2 - AL2  0.79862 0.147 417   5.422  <.0001
# AH2 - AL2  0.97554 0.329 296   2.968  0.0502
# AH3 - AL2  1.27689 0.323 321   3.950  0.0018

# Why not AA2 and AA1

lmer2.z <- lmer(z_scale ~ Taxa * Depth + (1|Loc/Depth), data = metadata)
summary(lmer2.z)
emmeans(lmer2.z, list(pairwise~Taxa:Depth), adjust = "tukey")
# AA1 20 - AA2 20            -0.4538 0.183 489.94  -2.475  0.6285 ??
# AA2 20 - AL2 20             0.8676 0.160 477.01   5.407  <.0001

k_boxplot(metadata, metadata$Taxa, metadata$z_scale, "Local height", colours, facet = "Depth")
# a bit different to outcrop proportions
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "Depth")


## John's comments ####
# sound models for design: m2 <- lmer(niche_variable ~ Taxa + (1|Loc/Depth), data = metadata)
# significant effect of taxa indicate differances in the average value of "niche variable" accounting for
# clustering within location and depth transect therein.
# To see if the differences in the niche variable depth on depth
# m_depth <- lmer(niche_variable ~ Taxa + Depth + Taxa:Depth + (1|Loc/Depth), data = metadata)

# multi-dimensional niche  ####
## create environmental distance matrix
AA_data <- metadata[metadata$Taxa == "AA1" | metadata$Taxa == "AA2",]
AA_env_dist_mat <- dist(AA_data[,c(12, 19:22)])

pcoa.AA <- pcoa(AA_env_dist_mat) # cmdscale(a.dist) same thing
pcoa.AA.dat <- as.data.frame(pcoa.AA$vectors[,1:5])
pcoa.AA.dat <- cbind(pcoa.AA.dat, AA_data$Taxa, AA_data$Depth)
ggplot(pcoa.AA.dat, aes(Axis.1, Axis.2,
                        color = `AA_data$Depth`, shape = `AA_data$Taxa`)) + geom_point()

# AA_env <- AA_data[,c(12,19:22)]
## fit model
AA_adonis <- adonis2(formula = AA_env_dist_mat ~ AA_data$Taxa, permutations = 999, strata = AA_data$Loc)
AA_adonis

AA_adonis_depth <- adonis2(formula = AA_env_dist_mat ~ AA_data$Taxa * AA_data$Depth, permutations = 999, strata = AA_data$Loc)
AA_adonis_depth
#AA_data$Depth                2     6779 0.03225 6.2874  0.003 **
#AA_data$Taxa:AA_data$Depth   1     1961 0.00933 3.6380  0.021 * 
# Something going on!!!


nmds <- metaMDS(comm = AA_data[,c(12, 19:22)], distance = "euclidean")
# Stress is nearly 0
library(viridis)

# Now we can build the plot!
dat <- data.frame(nmds$points)
dat$Taxa <- AA_data$Taxa
dat$Depth <- AA_data$Depth
ggplot() +
  geom_point(data = dat, aes(x = MDS1, y = MDS2,
                                     color = Taxa, shape = Depth), size = 3) +
  scale_color_manual(values = colours[1:2],
                     name = "Taxa") +
  annotate(geom = "label", x = -1, y = 1.25, size = 10,
           label = paste("Stress: ", round(nmds$stress, digits = 3))) +
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(size = 24))

## lamarcki ####
AL_data <- metadata[metadata$Taxa == "AL1" | metadata$Taxa == "AL2",]
AL_env_dist_mat <- dist(AL_data[,c(12, 19:22)])

AL_adonis <- adonis2(formula = AL_env_dist_mat ~ AL_data$Taxa, permutations = 999, strata = AL_data$Loc)
AL_adonis

AL_adonis_depth <- adonis2(formula = AL_env_dist_mat ~ AL_data$Taxa * AL_data$Depth, permutations = 999, strata = AL_data$Loc)
AL_adonis_depth

AH_data <- metadata[metadata$Taxa == "AH1" | metadata$Taxa == "AH2" | metadata$Taxa == "AH3",]
AH_env_dist_mat <- dist(AH_data[,c(12, 19:22)])

AH_adonis_depth <- adonis2(formula = AH_env_dist_mat ~ AH_data$Taxa * AH_data$Depth, permutations = 999, strata = AH_data$Loc)
AH_adonis_depth

# Just comparing depth
AA_data_d <- AA_data[AA_data$Depth == "20",]
AA_env_dist_mat_d <- dist(AA_data_d[,c(12, 19:22)])
AA_adonis_d <- adonis2(formula = AA_env_dist_mat_d ~ AA_data_d$Taxa, permutations = 999, strata = AA_data_d$Loc)
AA_adonis_d


# Single taxon ####
## AA2
k_boxplot(AA_data[AA_data$Taxa == "AA2",], AA_data$Depth[AA_data$Taxa == "AA2"],
          AA_data$ground_elevation_corr[AA_data$Taxa == "AA2"], "Ground Elevation corr", colours[1:3])
k_boxplot(AA_data[AA_data$Taxa == "AA2",], AA_data$Depth[AA_data$Taxa == "AA2"],
          AA_data$overhang_prop_2[AA_data$Taxa == "AA2"], "Overhang proportion", colours[1:3])
k_boxplot(AA_data[AA_data$Taxa == "AA2",], AA_data$Depth[AA_data$Taxa == "AA2"],
          AA_data$outcrop_prop[AA_data$Taxa == "AA2"], "Outcrop proportion", colours[1:3])
k_boxplot(AA_data[AA_data$Taxa == "AA2",], AA_data$Depth[AA_data$Taxa == "AA2"],
          AA_data$environment_rugosity_sqrt[AA_data$Taxa == "AA2"], "Environment rugosity", colours[1:3])
k_boxplot(AA_data[AA_data$Taxa == "AA2",], AA_data$Depth[AA_data$Taxa == "AA2"],
          AA_data$z_scale[AA_data$Taxa == "AA2"], "Local height", colours[1:3])
k_boxplot(AA_data[AA_data$Taxa == "AA2",], AA_data$Depth[AA_data$Taxa == "AA2"],
          AA_data$colony_rugosity[AA_data$Taxa == "AA2"], "Colony rugosity", colours[1:3])
k_boxplot(AA_data[AA_data$Taxa == "AA1",], AA_data$Depth[AA_data$Taxa == "AA1"],
          AA_data$colony_rugosity[AA_data$Taxa == "AA1"], "Colony rugosity", colours[1:3])
## AH1
k_boxplot(AH_data[AH_data$Taxa == "AH1",], AH_data$Depth[AH_data$Taxa == "AH1"],
          AH_data$ground_elevation_corr[AH_data$Taxa == "AH1"], "Ground Elevation corr", colours[1:3])
k_boxplot(AH_data[AH_data$Taxa == "AH1",], AH_data$Depth[AH_data$Taxa == "AH1"],
          AH_data$overhang_prop_2[AH_data$Taxa == "AH1"], "Overhang proportion", colours[1:3])
k_boxplot(AH_data[AH_data$Taxa == "AH1",], AH_data$Depth[AH_data$Taxa == "AH1"],
          AH_data$outcrop_prop[AH_data$Taxa == "AH1"], "Outcrop proportion", colours[1:3])
k_boxplot(AH_data[AH_data$Taxa == "AH1",], AH_data$Depth[AH_data$Taxa == "AH1"],
          AH_data$environment_rugosity_sqrt[AH_data$Taxa == "AH1"], "Environment rugosity", colours[1:3])
k_boxplot(AH_data[AH_data$Taxa == "AH1",], AH_data$Depth[AH_data$Taxa == "AH1"],
          AH_data$z_scale[AH_data$Taxa == "AH1"], "Local height", colours[1:3])
