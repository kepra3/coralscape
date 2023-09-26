# Title: Assessing differentiation of microhabitat by taxa using linear models
# Author: Katharine Prata
# Date: 21/6/22
# Last edit: 18/4/23

# Packages ####
library(tidyverse)
library(lme4)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(vegan)
library(MASS)
library(FSA)
library(viridis)
library(report)
library(ggplot2)

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
    p <- ggplot(df, aes(x = categorical, y = continuous, fill = categorical)) +
      geom_boxplot(fill = NA) +
      geom_point(size = 2, shape = 21) +
      #scale_x_discrete(labels = label) +
      scale_fill_manual(name = "Taxa", values = colours,
                         label = present_taxa) +
      theme_classic(base_size = 10) +
      ylab(lab) +
      xlab('Taxa') + 
      facet_wrap(~ Depth + Loc, nrow = 1)
  } else if (facet == "Depth") {
    p <- ggplot(df, aes(x = categorical, y = continuous, fill = categorical)) +
      geom_boxplot(fill = NA) +
      geom_point(size = 2, shape = 21) +
      #scale_x_discrete(labels = label) +
      scale_fill_manual(name = "Taxa", values = colours,
                         label = present_taxa) +
      theme_classic(base_size = 10) +
      ylab(lab) +
      xlab('Taxa') + 
      facet_wrap(~ Depth)
  } else if (facet == "Loc") {
    p <- ggplot(df, aes(x = categorical, y = continuous, fill = categorical)) +
      geom_boxplot(fill = NA) +
      geom_point(size = 2, shape = 21) +
      #scale_x_discrete(labels = label) +
      scale_fill_manual(name = "Taxa", values = colours,
                         label = present_taxa) +
      theme_classic(base_size = 10) +
      ylab(lab) +
      xlab('Taxa') + 
      facet_wrap(~ Loc)
  } else {
    label = paste0(present_taxa, "\n n = ", sample.sizes)
    p <- ggplot(df, aes(x = categorical, y = continuous, fill = categorical)) +
      geom_boxplot(fill = NA) +
      geom_point(size = 2, shape = 21) +
      scale_x_discrete(labels = label) +
      scale_fill_manual(name = "Taxa", values = colours,
                         label = present_taxa) +
      theme_classic(base_size = 10) +
      ylab(lab) +
      xlab('Taxa') }
  return(p)
}
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
mds_plot <- function(dat, cent, taxa.colours, ef.data) {
  plot <- ggplot() +
    geom_point(data = dat, aes(x = MDS1, y = MDS2, fill = Taxa, shape = Taxa), size = 1.5, alpha = 0.5) +
    geom_point(data = cent, aes(NMDS1, NMDS2, fill = Taxa, shape = Taxa), size = 3) +
    geom_point(data = cent, aes(NMDS1, NMDS2), size = 3, shape = 8, colour = "black") +
    scale_fill_manual(values = taxa.colours, name = "Taxa") +
    scale_shape_manual(values = c(22, 24), name = "Taxa") +
    #coord_fixed() + ## need aspect ratio of 1
    ylim(c(-4.5, 4.5)) +
    xlim(c(-4.5, 4.5)) +
    #ggtitle(paste("Stress: ", round(nmds$stress, digits = 3))) +
    theme_minimal() +
    geom_segment(data = ef.data, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                 arrow = arrow(length = unit (0.25, "cm")), color = "blue") +
    geom_text(data = ef.data, aes(x = NMDS1, y = NMDS2, label = row.names(ef.data)),
              colour = "blue", size = 3) +
    theme(legend.position = "none",
          text = element_text(size = 8),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.ticks = element_line(),
          axis.text = element_text(size = 8))
  return(plot)
}

# Import datasets
setwd("~/git/coralscape/results/microhabitat_paper_results/")
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
metadata$Loc <- factor(metadata$Loc, levels = c("WP", "SB", "CA", "SQ"))
metadata$Site <- factor(metadata$Site, levels = c("WP05", "CA05", "SB05", "SQ05",
                                                  "WP10", "CA10", "SB10", "SQ12", 
                                                  "WP20", "CA20", "SB20", "SQ20"))

                       
metadata[, c(2, 4:6)] <- data.frame(lapply(metadata[, c(2, 4:6)], as.factor))
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
site_list <- levels(metadata$Site)
for  (i in site_list) {
  metadata$z_scale[metadata$Site == i] <- scale(metadata$z[metadata$Site == i])
  
}

# Separate by Depths ####
metadata20 <- metadata[metadata$Depth == "20",]
metadata10 <- metadata[metadata$Depth == "10-12" &
                         (metadata$Taxa == "AA2" |
                            metadata$Taxa == "AH1" |
                            metadata$Taxa == "AH3" |
                            metadata$Taxa == "AL2"),]
metadata5 <- metadata[metadata$Depth == "5",]


#### ********************** STOP HERE IF WANT TO SKIP TO ANOTHER SECTION ************************** ####

# Proportion of Taxa by Depth ####
library(corrplot)
tab <- table(metadata$Taxa, metadata$Depth)
# relativised by abudance due to bias in total taxa
tab.stand <- round((tab/rowSums(tab)) * 100, 0)
chi <- chisq.test(tab.stand)
fish <- fisher.test(table(metadata$Taxa, metadata$Depth), simulate.p.value = TRUE)
corrplot(chi$residuals, is.corr = FALSE, col.lim = c(-11, 11))
chi$expected
chi$observed

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

# Depth and location ####
x <- as.data.frame(table(metadata$Taxa, metadata$Site))
str(x)
x$Freq <- as.numeric(x$Freq)
x$bin[x$Freq == 0] <- NA
x$bin[x$Freq >= 3 & x$Freq < 10] <- 4
x$bin[x$Freq >= 10 & x$Freq < 20] <- 5
x$bin[x$Freq >= 20 & x$Freq < 30] <- 6
x$bin[x$Freq >= 30 & x$Freq < 40] <- 7
x$bin[x$Freq >= 40 & x$Freq < 50] <- 8
x$bin[x$Freq >= 50] <- 9
x$Freq[x$Freq == 0] <- NA
ggplot(data = x, aes(Var1, 1, fill = Var1)) + geom_point(size = x$bin, shape = 21) +
  scale_fill_manual(values = colours) + 
  geom_text(aes(Var1, 1, label = Freq), size = 2, colour = "black", face = "bold") +
  facet_wrap(~Var2, nrow = 4, ncol = 4) + theme_classic() +
  scale_x_discrete(position = "top") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        strip.placement = "outside",
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA),
        panel.grid.major.x = element_line(),
        legend.position = "none")
#ggsave("depth_distribution.pdf", height = 4, width = 17, units = "cm")





# *** Linear Models ***  ####
# Phenotypic differences ####
## Colony rugosity ####
k_boxplot(metadata5, metadata5$Taxa, metadata5$colony_rugosity, "Colony rugosity", colours[2:5], facet = "All")
k_boxplot(metadata10, metadata10$Taxa, metadata10$colony_rugosity, "Colony rugosity", colours[c(2,3,5,7)], facet = "All")
k_boxplot(metadata20, metadata20$Taxa, metadata20$colony_rugosity, "Colony rugosity", colours[c(1:3,6,7)], facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$colony_rugosity, "Colony rugosity", colours, facet = "Depth")
# ggsave("ground_elevation_by_site.pdf", height = 11, width = 20, units = "cm")

# Random effect of Location
lmer.colony_rugosity20 <- lmer(colony_rugosity ~ Taxa + (1|Loc), data = metadata20)
summary(lmer.colony_rugosity20)
plot(lmer.colony_rugosity20)
emmeans(lmer.colony_rugosity20, list(pairwise ~ Taxa), adjust = "tukey")
# AA1 - AL1   0.2964 0.0858 300   3.452  0.0057
# AA1 - AL2   0.2638 0.0734 299   3.596  0.0034
# AA2 - AL1   0.2652 0.0687 299   3.860  0.0013
# AA2 - AL2   0.2326 0.0520 300   4.472  0.0001
# conclusion: AA1/AA2 are more rugose than AL1/AL2

# Random effect of Location
lmer.colony_rugosity10 <- lmer(colony_rugosity ~ Taxa + (1|Loc), data = metadata10)
summary(lmer.colony_rugosity10)
plot(lmer.colony_rugosity10)
emmeans(lmer.colony_rugosity10, list(pairwise ~ Taxa), adjust = "tukey")

# Random effect of Location
lmer.colony_rugosity5 <- lmer(colony_rugosity ~ Taxa + (1|Loc), data = metadata5)
summary(lmer.colony_rugosity5)
plot(lmer.colony_rugosity5)
emmeans(lmer.colony_rugosity5, list(pairwise ~ Taxa), adjust = "tukey")

## Colony radius####
k_boxplot(metadata5, metadata5$Taxa, metadata5$colony_range, "Colony radius", colours[2:5], facet = "All")
k_boxplot(metadata10, metadata10$Taxa, metadata10$colony_range, "Colony radius", colours[c(2,3,5,7)], facet = "All")
k_boxplot(metadata20, metadata20$Taxa, metadata20$colony_range, "Colony radius", colours[c(1:3,6,7)], facet = "All")
k_boxplot(metadata, metadata$Taxa, metadata$colony_range, "Colony radius", colours, facet = "Depth")
# ggsave("ground_elevation_by_site.pdf", height = 11, width = 20, units = "cm")

# Random effect of Location
lmer.colony_radius20 <- lmer(colony_range ~ Taxa + (1|Loc), data = metadata20)
summary(lmer.colony_radius20)
plot(lmer.colony_radius20)
emmeans(lmer.colony_radius20, list(pairwise ~ Taxa), adjust = "tukey")
# AA1 - AL1 -0.08907 0.01387 299  -6.421  <.0001
# AA1 - AL2 -0.06300 0.01186 300  -5.312  <.0001
# AA2 - AL1 -0.08219 0.01110 298  -7.402  <.0001
# AA2 - AL2 -0.05612 0.00841 300  -6.673  <.0001
# AH1 - AL1 -0.11977 0.02184 299  -5.484  <.0001
# AH1 - AL2 -0.09370 0.02057 299  -4.556  0.0001
# conclusion: AA1/AA2 are more rugose than AL1/AL2

# Random effect of Location
lmer.colony_radius10 <- lmer(colony_range ~ Taxa + (1|Loc), data = metadata10)
summary(lmer.colony_radius10)
plot(lmer.colony_radius10)
emmeans(lmer.colony_radius10, list(pairwise ~ Taxa), adjust = "tukey")
# AA2 - AH1   0.0562 0.0141 148   3.997  0.0006
# AA2 - AH3   0.0698 0.0229 148   3.045  0.0145
# AH1 - AL2  -0.0869 0.0227 148  -3.832  0.0011
# AH3 - AL2  -0.1005 0.0290 148  -3.462  0.0039

# Random effect of Location
lmer.colony_radius5 <- lmer(colony_range ~ Taxa + (1|Loc), data = metadata5)
summary(lmer.colony_radius)
plot(lmer.colony_radius5)
emmeans(lmer.colony_radius5, list(pairwise ~ Taxa), adjust = "tukey")

# are colony rugosity and range correlated

lmer.colony_range.radius <- lmer(colony_range ~ colony_rugosity + (1|Taxa), data = metadata5)
summary(lmer.colony_radius.radius)
emmeans(lmer.colony_radius/radius, ~ colony_rugosity)
# ****** Microhabitat differences ******####
## Ground elevation ####

a5 <- k_boxplot(metadata5, metadata5$Taxa, metadata5$ground_elevation_corr, "Substrate angle", colours[2:5], facet = "All")
ggsave("new_plots/ground_elevation_corr_5m.pdf", height = 10, width = 20, units = "cm")
a10 <- k_boxplot(metadata10, metadata10$Taxa, metadata10$ground_elevation_corr, "Substrate angle", colours[c(2,3,5,7)], facet = "All")
ggsave("new_plots/ground_elevation_corr_10m.pdf", height = 10, width = 20, units = "cm")
a20 <- k_boxplot(metadata20, metadata20$Taxa, metadata20$ground_elevation_corr, "Substrate angle", colours[c(1:3,6,7)], facet = "All")
ggsave("new_plots/ground_elevation_corr_20m.pdf", height = 10, width = 20, units = "cm")
k_boxplot(metadata, metadata$Taxa, metadata$ground_elevation_corr, "Substrate angle", colours, facet = "Depth")
ggsave("new_plots/ground_elevation_corr_all.pdf", height = 10, width = 20, units = "cm")


figure <- ggarrange(a5, a10, a20,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
figure
ggsave("new_plots/ch4_S1_substrate_angle.pdf", height = 26, width = 15, units = "cm")
# ggsave("ground_elevation_by_site.pdf", height = 11, width = 20, units = "cm")

# Random effect of Location
lmer.ground_elevation20 <- lmer(ground_elevation_corr ~ Taxa + (1|Loc), data = metadata20)
summary(lmer.ground_elevation20)
plot(lmer.ground_elevation20)
emmeans(lmer.ground_elevation20, list(pairwise ~ Taxa), adjust = "tukey")
# AA2 - AL2    11.64 3.71 299   3.141  0.0158
ground_elevation20.emmeans <- emmeans(lmer.ground_elevation20, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(ground_elevation20.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/ground_elevation20.emmeans.txt", quote = FALSE)

# conclusion: AA2 sits at higher angles than AL2 at 20m

# Random effect of Location
lmer.ground_elevation10 <- lmer(ground_elevation_corr ~ Taxa + (1|Loc), data = metadata10)
summary(lmer.ground_elevation10)
plot(lmer.ground_elevation10)
emmeans(lmer.ground_elevation10, list(pairwise ~ Taxa), adjust = "tukey")
ground_elevation10.emmeans <- emmeans(lmer.ground_elevation10, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(ground_elevation10.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/ground_elevation10.emmeans.txt", quote = FALSE)

# Random effect of Location
lmer.ground_elevation5 <- lmer(ground_elevation_corr ~ Taxa + (1|Loc), data = metadata5)
summary(lmer.ground_elevation5)
plot(lmer.ground_elevation5)
emmeans(lmer.ground_elevation5, list(pairwise ~ Taxa), adjust = "tukey")
ground_elevation5.emmeans <- emmeans(lmer.ground_elevation5, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(ground_elevation5.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/ground_elevation5.emmeans.txt", quote = FALSE)

### is colony rugosity affecting ground elevation ####
ggplot(metadata20, aes(ground_elevation_corr, colony_rugosity, colour = Taxa)) + geom_point() + geom_smooth(method = "lm")
ggplot(metadata20, aes(ground_elevation_corr, environment_rugosity_sqrt, colour = Taxa)) + geom_point() + geom_smooth(method = "lm")
lmer.rugosity_elevation20 <- lmer(ground_elevation_corr ~ colony_rugosity + (1|Taxa), data = metadata20)
summary(lmer.rugosity_elevation20)
plot(lmer.rugosity_elevation20)
anova(lmer.rugosity_elevation20)

## Environment rugosity ####
r5 <- k_boxplot(metadata5, metadata5$Taxa, metadata5$environment_rugosity_sqrt, "Square-root of environment rugosity", colours[2:5], facet = "All")
ggsave("new_plots/environment_rugosity_sqrt_5m.pdf", height = 10, width = 20, units = "cm")
r10 <- k_boxplot(metadata10, metadata10$Taxa, metadata10$environment_rugosity_sqrt, "Square-root of environment rugosity", colours[c(2,3,5,7)], facet = "All")
ggsave("new_plots/environment_rugosity_sqrt_10m.pdf", height = 10, width = 20, units = "cm")
r20 <- k_boxplot(metadata20, metadata20$Taxa, metadata20$environment_rugosity_sqrt, "Square-root of environment rugosity", colours[c(1:3,6,7)], facet = "All")
ggsave("new_plots/environment_rugosity_sqrt_20m.pdf", height = 10, width = 20, units = "cm")
k_boxplot(metadata, metadata$Taxa, metadata$environment_rugosity_sqrt, "Square-root of environment rugosity", colours, facet = "Depth")
ggsave("new_plots/environment_rugosity_sqrt_all.pdf", height = 10, width = 20, units = "cm")

# supp figure
figure <- ggarrange(r5, r10, r20,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
figure
ggsave("new_plots/ch4_S1_environment_rugosity.pdf", height = 26, width = 15, units = "cm")


# Including Loc/Depth as random effect
lmer.environment_rugosity5 <- lmer(environment_rugosity_sqrt ~ Taxa + (1|Loc), data = metadata5)
summary(lmer.environment_rugosity5)
plot(lmer.environment_rugosity5)
emmeans(lmer.environment_rugosity5, list(pairwise ~ Taxa), adjust = "tukey")
environment_rugosity5.emmeans <- emmeans(lmer.environment_rugosity5, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(environment_rugosity5.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/environment_rugosity5.emmeans.txt", quote = FALSE)

# Including Loc/Depth as random effect
lmer.environment_rugosity10 <- lmer(environment_rugosity_sqrt ~ Taxa + (1|Loc), data = metadata10)
summary(lmer.environment_rugosity10)
plot(lmer.environment_rugosity10)
emmeans(lmer.environment_rugosity10, list(pairwise ~ Taxa), adjust = "tukey")
# AA2 - AL2   0.1511 0.0595 148   2.540  0.0580
# AH1 - AL2   0.2009 0.0725 148   2.771  0.0316 *
# AH3 - AL2   0.2252 0.0928 148   2.427  0.0765
environment_rugosity10.emmeans <- emmeans(lmer.environment_rugosity10, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(environment_rugosity10.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/environment_rugosity10.emmeans.txt", quote = FALSE)

# Including Loc/Depth as random effect
lmer.environment_rugosity20 <- lmer(environment_rugosity_sqrt ~ Taxa + (1|Loc), data = metadata20)
summary(lmer.environment_rugosity20)
plot(lmer.environment_rugosity20)
emmeans(lmer.environment_rugosity20, list(pairwise ~ Taxa), adjust = "tukey")
# AA2 - AL1  0.09492 0.0302 297   3.140  0.0159
# AA2 - AL2  0.06703 0.0229 298   2.928  0.0299
environment_rugosity20.emmeans <- emmeans(lmer.environment_rugosity20, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(environment_rugosity20.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/environment_rugosity20.emmeans.txt", quote = FALSE)

# Conclusions: 

## Outcrop proportion ####
out5 <- k_boxplot(metadata5, metadata5$Taxa, metadata5$outcrop_prop, "Outcrop proportion", colours[2:5], facet = "All")
ggsave("new_plots/outcrop_prop_5m.pdf", height = 10, width = 20, units = "cm")
out10 <- k_boxplot(metadata10, metadata10$Taxa, metadata10$outcrop_prop, "Outcrop proportion", colours[c(2,3,5,7)], facet = "All")
ggsave("new_plots/outcrop_prop_10m.pdf", height = 10, width = 20, units = "cm")
out20 <- k_boxplot(metadata20, metadata20$Taxa, metadata20$outcrop_prop, "Outcrop proportion", colours[c(1:3,6,7)], facet = "All")
ggsave("new_plots/outcrop_prop_20m.pdf", height = 10, width = 20, units = "cm")
k_boxplot(metadata, metadata$Taxa, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "Depth")
ggsave("new_plots/outcrop_prop_all.pdf", height = 10, width = 20, units = "cm")

# supp figure
figure <- ggarrange(out5, out10, out20,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
figure
ggsave("new_plots/ch4_S3_outcrop_proportion.pdf", height = 26, width = 15, units = "cm")



# Including Loc/Depth as random effect
lmer.outcrop_prop5 <- lmer(outcrop_prop ~ Taxa + (1|Loc), data = metadata5)
summary(lmer.outcrop_prop5)
emmeans(lmer.outcrop_prop5, list(pairwise ~ Taxa), adjust = "tukey")
outcrop_prop5.emmeans <- emmeans(lmer.outcrop_prop5, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(outcrop_prop5.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/outcrop_prop5.emmeans.txt", quote = FALSE)

# Including Loc/Depth as random effect
lmer.outcrop_prop10 <- lmer(outcrop_prop ~ Taxa + (1|Loc), data = metadata10)
summary(lmer.outcrop_prop10)
emmeans(lmer.outcrop_prop10, list(pairwise ~ Taxa), adjust = "tukey")
outcrop_prop10.emmeans <- emmeans(lmer.outcrop_prop10, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(outcrop_prop10.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/outcrop_prop10.emmeans.txt", quote = FALSE)


# Including Loc/Depth as random effect
lmer.outcrop_prop20 <- lmer(outcrop_prop ~ Taxa + (1|Loc), data = metadata20)
summary(lmer.outcrop_prop20)
emmeans(lmer.outcrop_prop20, list(pairwise ~ Taxa), adjust = "tukey")
#AA1 - AL2  0.12005 0.0277 299   4.329  0.0002
#AA2 - AL1  0.10298 0.0260 299   3.966  0.0009
#AA2 - AL2  0.15465 0.0197 300   7.866  <.0001
outcrop_prop20.emmeans <- emmeans(lmer.outcrop_prop20, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(outcrop_prop20.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/outcrop_prop20.emmeans.txt", quote = FALSE)

## Overhang proportion ####
# TODO: USE BETA REGRESSION THAT ACCOUNTS FOR ZERO INFLATED DATA ##
ov5 <- k_boxplot(metadata5, metadata5$Taxa, metadata5$overhang_prop_2, "Overhang proportion", colours[2:5], facet = "All")
ggsave("new_plots/overhang_prop_2_5m.pdf", height = 10, width = 20, units = "cm")
ov10 <- k_boxplot(metadata10, metadata10$Taxa, metadata10$overhang_prop_2, "Overhang proportion", colours[c(2,3,5,7)], facet = "All")
ggsave("new_plots/overhang_prop_2_10m.pdf", height = 10, width = 20, units = "cm")
ov20 <- k_boxplot(metadata20, metadata20$Taxa, metadata20$overhang_prop_2, "Overhang proportion", colours[c(1:3,6,7)], facet = "All")
ggsave("new_plots/overhang_prop_2_20m.pdf", height = 10, width = 20, units = "cm")
k_boxplot(metadata, metadata$Taxa, metadata$overhang_prop_2, "Overhang proportion", colours, facet = "Depth")
ggsave("new_plots/overhang_prop_2_all.pdf", height = 10, width = 20, units = "cm")

# supp figure
figure <- ggarrange(ov5, ov10, ov20,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
figure
ggsave("new_plots/ch4_S2_overhang_proportion.pdf", height = 26, width = 15, units = "cm")


# Including Loc/Depth as random effect
lmer.overhang_prop5 <- glmer(overhang_prop_2 ~ Taxa + (1|Loc), data = metadata5, family = binomial(link = "logit"))
summary(lmer.overhang_prop5)
plot(lmer.overhang_prop5)
overdisp_fun(lmer.overhang_prop5)
fm1.all <- allFit(lmer.overhang_prop5)
ss <- summary(fm1.all)
ss$ fixef               ## fixed effects
ss$ llik                ## log-likelihoods
ss$ sdcor               ## SDs and correlations
ss$ theta               ## Cholesky factors
ss$ which.OK 

emmeans(lmer.overhang_prop5, list(pairwise ~ Taxa), adjust = "tukey")
overhang_prop5.emmeans <- emmeans(lmer.overhang_prop5, list(pairwise ~ Taxa), adjust = "tukey")
plot(overhang_prop5.emmeans)



write.csv(overhang_prop5.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/overhang_prop5.emmeans.txt", quote = FALSE)

# Including Loc/Depth as random effect
lmer.overhang_prop10 <- glmer(overhang_prop_2 ~ Taxa + (1|Loc), data = metadata10, family = binomial(link = "logit"))
summary(lmer.overhang_prop10)
plot(lmer.overhang_prop10)
emmeans(lmer.overhang_prop10, list(pairwise ~ Taxa), adjust = "tukey")
overhang_prop10.emmeans <- emmeans(lmer.overhang_prop10, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(overhang_prop10.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/overhang_prop10.emmeans.txt", quote = FALSE)


# Including Loc/Depth as random effect
lmer.overhang_prop20 <- glmer(overhang_prop_2 ~ Taxa + (1|Loc), data = metadata20, family = binomial(link = "logit"))
summary(lmer.overhang_prop20)
plot(lmer.overhang_prop20)
emmeans(lmer.overhang_prop20, list(pairwise ~ Taxa), adjust = "tukey")
# AA2 - AL2  -1.4900 0.431 Inf  -3.457  0.0050
overhang_prop20.emmeans <- emmeans(lmer.overhang_prop20, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(overhang_prop20.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/overhang_prop20.emmeans.txt", quote = FALSE)

## Z - local depth... ####
z5 <- k_boxplot(metadata5, metadata5$Taxa, metadata5$z_scale, "Local height", colours[2:5], facet = "All")
ggsave("new_plots/z_scale_5m.pdf", height = 10, width = 20, units = "cm")
z10 <- k_boxplot(metadata10, metadata10$Taxa, metadata10$z_scale, "Local height", colours[c(2,3,5,7)], facet = "All")
ggsave("new_plots/z_scale_10m.pdf", height = 10, width = 20, units = "cm")
z20 <- k_boxplot(metadata20, metadata20$Taxa, metadata20$z_scale, "Local height", colours[c(1:3,6,7)], facet = "All")
ggsave("new_plots/z_scale_20m.pdf", height = 10, width = 20, units = "cm")
k_boxplot(metadata, metadata$Taxa, metadata$z_scale, "Local height", colours, facet = "Depth")
ggsave("new_plots/z_scale_all.pdf", height = 10, width = 20, units = "cm")

# supp figure
figure <- ggarrange(z5, z10, z20,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
figure
ggsave("new_plots/ch4_S4_local_height.pdf", height = 26, width = 15, units = "cm")


# Including Loc/Depth as random effect
lmer.z_scale5 <- lmer(z_scale ~ Taxa + (1|Loc), data = metadata5)
summary(lmer.z_scale5)
plot(lmer.z_scale5)
emmeans(lmer.z_scale5, list(pairwise ~ Taxa), adjust = "tukey")
z_scale5.emmeans <- emmeans(lmer.z_scale5, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(z_scale5.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/z_scale5.emmeans.txt", quote = FALSE)

# Including Loc/Depth as random effect
lmer.z_scale10 <- lmer(z_scale ~ Taxa + (1|Loc), data = metadata10)
summary(lmer.z_scale10)
plot(lmer.z_scale10)
emmeans(lmer.z_scale10, list(pairwise ~ Taxa), adjust = "tukey")
# AH3 - AL2    1.765 0.645 143   2.736  0.0350
z_scale10.emmeans <- emmeans(lmer.z_scale10, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(z_scale10.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/z_scale10.emmeans.txt", quote = FALSE)

# Including Loc/Depth as random effect
lmer.z_scale20 <- lmer(z_scale ~ Taxa + (1|Loc), data = metadata20)
summary(lmer.z_scale20)
plot(lmer.z_scale20)
emmeans(lmer.z_scale20, list(pairwise ~ Taxa), adjust = "tukey")
# AA1 - AA2  -0.4538 0.182 300  -2.491  0.0953
# AA2 - AL2   0.8676 0.159 294   5.455  <.0001
z_scale20.emmeans <- emmeans(lmer.z_scale20, list(pairwise ~ Taxa), adjust = "tukey")
write.csv(z_scale20.emmeans$`pairwise differences of Taxa`, file = "univariate_emmeans/z_scale20.emmeans.txt", quote = FALSE)

# Multi-dimensional niche  ####
## agaricites ####
## create environmental distance matrix
AA_data <- metadata20[metadata20$Taxa == "AA1" | metadata20$Taxa == "AA2",]
env <- AA_data[,c(12, 19, 20, 21, 22)]
env$outcrop_prop <- scale(env$outcrop_prop) 
env$outcrop_prop <- env$outcrop_prop + abs(min(env$outcrop_prop))
env$environment_rugosity_sqrt <- scale(env$environment_rugosity_sqrt)
env$environment_rugosity_sqrt <- env$environment_rugosity_sqrt + abs(min(env$environment_rugosity_sqrt))
env$ground_elevation_corr <- scale(env$ground_elevation_corr)
env$ground_elevation_corr <- env$ground_elevation_corr + abs(min(env$ground_elevation_corr))
env$overhang_prop_2 <- scale(env$overhang_prop_2)
env$overhang_prop_2 <-  env$overhang_prop_2+ abs(min(env$overhang_prop_2))
env$z_scale <- env$z_scale + abs(min(env$z_scale))
AA_env_dist_mat <- dist(env)

## fit model
AA_adonis <- adonis2(formula = AA_env_dist_mat ~ AA_data$Taxa, permutations = 999, strata = AA_data$Loc, by = "term")
AA_adonis
names(AA_adonis)

# Plot
nmds <- metaMDS(comm = env, distance = "euclidean", k = 2)
ef <- envfit(nmds, env, permu = 999, strata = AA_data$Loc)
ef.data <- as.data.frame(scores(ef, display = "vectors"))
row.names(ef.data) <- c("outcrop", "rugosity", "elevation", "overhang", "z")
dat <- data.frame(nmds$points)
dat$Taxa <- AA_data$Taxa
# Calculate centroids
scores <- scores(nmds, display = "sites")
cent <- aggregate(scores ~ Taxa, data = dat, FUN = "mean")
# Make plot
ac.mds <- mds_plot(dat, cent, colours[1:2], ef.data)
ac.mds
ggsave(plot = ac.mds, "NMDS_Agaricites.pdf", height = 8, width = 8, units = "cm")

## humilis ####
AH_data <- metadata5[metadata5$Taxa == "AH1" | metadata5$Taxa == "AH2" | metadata5$Taxa == "AH3",]
env <- AH_data[,c(12, 19, 20, 21, 22)]
env$outcrop_prop <- scale(env$outcrop_prop)
env$environment_rugosity_sqrt <- scale(env$environment_rugosity_sqrt)
env$ground_elevation_corr <- scale(env$ground_elevation_corr)
env$overhang_prop_2 <- scale(env$overhang_prop_2)
AH_env_dist_mat <- dist(env)
AH_adonis <- adonis2(formula = AH_env_dist_mat ~ AH_data$Taxa, permutations = 999, strata = AH_data$Loc)
AH_adonis

AH_data <- metadata5[metadata5$Taxa == "AH1" | metadata5$Taxa == "AH2",]
env <- AH_data[,c(12, 19, 20, 21, 22)]
env$outcrop_prop <- scale(env$outcrop_prop)
env$environment_rugosity_sqrt <- scale(env$environment_rugosity_sqrt)
env$ground_elevation_corr <- scale(env$ground_elevation_corr)
env$overhang_prop_2 <- scale(env$overhang_prop_2)
AH_env_dist_mat <- dist(env)
AH_adonis <- adonis2(formula = AH_env_dist_mat ~ AH_data$Taxa, permutations = 999, strata = AH_data$Loc)
AH_adonis

AH_data <- metadata5[metadata5$Taxa == "AH1" | metadata5$Taxa == "AH3",]
env <- AH_data[,c(12, 19, 20, 21, 22)]
env$outcrop_prop <- scale(env$outcrop_prop)
env$environment_rugosity_sqrt <- scale(env$environment_rugosity_sqrt)
env$ground_elevation_corr <- scale(env$ground_elevation_corr)
env$overhang_prop_2 <- scale(env$overhang_prop_2)
AH_env_dist_mat <- dist(env)
AH_adonis <- adonis2(formula = AH_env_dist_mat ~ AH_data$Taxa, permutations = 999, strata = AH_data$Loc)
AH_adonis

AH_data <- metadata5[metadata5$Taxa == "AH2" | metadata5$Taxa == "AH3",]
env <- AH_data[,c(12, 19, 20, 21, 22)]
env$outcrop_prop <- scale(env$outcrop_prop)
env$environment_rugosity_sqrt <- scale(env$environment_rugosity_sqrt)
env$ground_elevation_corr <- scale(env$ground_elevation_corr)
env$overhang_prop_2 <- scale(env$overhang_prop_2)
AH_env_dist_mat <- dist(env)
AH_adonis <- adonis2(formula = AH_env_dist_mat ~ AH_data$Taxa, permutations = 999, strata = AH_data$Loc)
AH_adonis

## lamarcki ####
AL_data <- metadata20[metadata20$Taxa == "AL1" | metadata20$Taxa == "AL2",]
env <- AL_data[,c(12, 19, 20, 21, 22)]
env$outcrop_prop <- scale(env$outcrop_prop) 
env$outcrop_prop <- env$outcrop_prop + abs(min(env$outcrop_prop))
env$environment_rugosity_sqrt <- scale(env$environment_rugosity_sqrt)
env$environment_rugosity_sqrt <- env$environment_rugosity_sqrt + abs(min(env$environment_rugosity_sqrt))
env$ground_elevation_corr <- scale(env$ground_elevation_corr)
env$ground_elevation_corr <- env$ground_elevation_corr + abs(min(env$ground_elevation_corr))
env$overhang_prop_2 <- scale(env$overhang_prop_2)
env$overhang_prop_2 <-  env$overhang_prop_2+ abs(min(env$overhang_prop_2))
env$z_scale <- env$z_scale + abs(min(env$z_scale))
AL_env_dist_mat <- dist(env)
AL_adonis <- adonis2(formula = AL_env_dist_mat ~ AL_data$Taxa, permutations = 999, strata = AL_data$Loc)


AL_adonis <- adonis2(formula = AL_env_dist_mat ~ AL_data$Taxa,
                     permutations = 999, strata = AL_data$Loc, by = "terms")
AL_adonis
# the adonis test are identical to the annova.cca of dbrda
# permutational MANOVA
nmds <- metaMDS(env, k=2, distance = "euclidean") # K maybe needs to be 3
ef <- envfit(nmds, env, permu = 999, strata = AL_data$Loc)

ef.data <- as.data.frame(scores(ef, display = "vectors"))
row.names(ef.data) <- c("outcrop", "rugosity", "elevation", "overhang", "z")
dat <- data.frame(nmds$points)
dat$Taxa <- AL_data$Taxa
# Calculate centroids
scores <- scores(nmds, display = "sites")
cent <- aggregate(scores ~ Taxa, data = dat, FUN = "mean")

lm.mds <- mds_plot(dat, cent, colours[6:7], ef.data)
lm.mds
# Plot
ggsave(plot = lm.mds, "NMDS_Lamarcki.pdf", height = 8, width = 8, units = "cm")

