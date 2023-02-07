# Packages
library(tidyverse)
library(lme4)
library(emmeans)
library(ggplot2)
library(vegan)
library(MASS)

# Functions
combine_metadata <- function(struc.complex, clusters, coordinates) {
  # match
  struc.complex.subset <- subset(struc.complex, sample_name %in% clusters$Individual)
  clusters.subset <- subset(clusters, Individual %in% struc.complex.subset$sample_name)
  coordinates.subset <- subset(coordinates, X %in% struc.complex.subset$sample_name)
  
  # sort
  sorted.clusters.subset <- clusters.subset[order(clusters.subset[,1]),]
  sorted.coordinates.subset <- coordinates.subset[order(coordinates.subset[,1]),]
  sorted.struc.complex.subset <- struc.complex.subset[order(struc.complex.subset[,2]),]
  # check if identiical
  cat("\nChecking clusters and struc complex: ", identical(sorted.clusters.subset[,1], sorted.struc.complex.subset[,2]))
  cat("\nChecking clusters and coordinates: ", identical(sorted.clusters.subset[,1], sorted.coordinates.subset[,1]))
  # Put metadata together
  metadata <- cbind(sorted.clusters.subset, sorted.struc.complex.subset[,3:10], sorted.coordinates.subset[,-1])
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
  present_clusters <- levels(categorical)
  sample.sizes <- NA
  
  for (i in 1:length(present_clusters)) {
    sample.sizes[i] <- length(categorical[categorical == present_clusters[i]])
  }
  
  
  if (facet == "All") {
  p <- ggplot(df, aes(x = categorical, y = continuous, color = categorical)) +
    geom_boxplot(notch = TRUE) +
    geom_point(size = 3, shape = 1) +
    #scale_x_discrete(labels = label) +
    scale_color_manual(name = "Clusters", values = colours,
                       label = present_clusters) +
    theme_classic(base_size = 12) +
    ylab(lab) +
    xlab('Clusters') + 
    facet_wrap(~ Depth + Loc)
  } else if (facet == "Depth") {
    p <- ggplot(df, aes(x = categorical, y = continuous, color = categorical)) +
    geom_boxplot(notch = TRUE) +
    geom_point(size = 3, shape = 1) +
    #scale_x_discrete(labels = label) +
    scale_color_manual(name = "Clusters", values = colours,
                       label = present_clusters) +
    theme_classic(base_size = 12) +
    ylab(lab) +
    xlab('Clusters') + 
    facet_wrap(~ Depth)
  } else if (facet == "Loc") {
    p <- ggplot(df, aes(x = categorical, y = continuous, color = categorical)) +
      geom_boxplot(notch = TRUE) +
      geom_point(size = 3, shape = 1) +
      #scale_x_discrete(labels = label) +
      scale_color_manual(name = "Clusters", values = colours,
                         label = present_clusters) +
      theme_classic(base_size = 12) +
      ylab(lab) +
      xlab('Clusters') + 
      facet_wrap(~ Loc)
  } else {
    label = paste0(present_clusters, "\n n = ", sample.sizes)
    p <- ggplot(df, aes(x = categorical, y = continuous, colour = categorical)) +
    geom_boxplot(notch = TRUE) +
    geom_point(size = 3, shape = 1) +
    scale_x_discrete(labels = label) +
    scale_color_manual(name = "Clusters", values = colours,
                         label = present_clusters) +
    theme_classic(base_size = 12) +
    ylab(lab) +
    xlab('Clusters') }
  return(p)
}
k_rda <- function(env.scale, depth) {
  env.scale$Site <- metadata$Site
  env.scale$Depth <- metadata$Depth
  env.scale.depth <- env.scale[env.scale$Depth == depth,]
  # see relationships
  pairs(env.scale.depth, lower.panel = NULL, col = as.numeric(env.scale.depth$Clusters))
  rda <- rda(env.scale.depth[,-c(5,6,7)], scale = TRUE)
  #make plot
  biplot(rda,
         display = c("sites", 
                     "species"),
         type = c("text",
                  "points"))
  #Add "hulls"
  env.scale.depth$Clusters <- droplevels(env.scale.depth$Clusters)
  cluster.depth.names <- levels(env.scale.depth$Clusters)
  ordihull(rda,
           group = env.scale.depth$Clusters, col = rep(1:length(cluster.depth.names)))
  legend("topleft",
         col = rep(1:length(cluster.depth.names)),
         lty = 1,
         legend = cluster.depth.names)
  title(paste("PCA", depth, "m"))
  colvec <- rep(1:length(cluster.depth.names))
  with(env.scale.depth, points(rda, display = "sites", col = colvec[Clusters], pch = 21, bg = colvec[Clusters]))
  #scores <- data.frame(rda$CA$u)
  #scores$group <- env.scale.depth$Clusters
  #pca.centroids <- aggregate(scores[,1:4], list(scores[,5]), mean)
  #points(pca.centroids$PC1*4, pca.centroids$PC2*4, col = colvec)
  print('Displaying plot')
}
k_lda <- function(env.scale.lda, species) {
  env.scale.species <- env.scale.lda[env.scale.lda$Clusters == paste0(species, 1) | env.scale.lda$Clusters == paste0(species, 2),]
  env.scale.species$Clusters <- droplevels(env.scale.species$Clusters)
  env.scale.species.lda <- env.scale.species[,-7]
  lda <- lda(formula = Clusters ~ ., data = env.scale.species.lda)
  
  ## get the x,y coordinates for the LDA plot
  data.lda.values <- predict(lda)
  
  ## create a dataframe that has all the info we need to draw a graph
  plot.data <- data.frame(X = data.lda.values$x[,1], Clusters = env.scale.species$Clusters, Depth = env.scale.species$Depth)
  
  head(plot.data)
  
  ## draw a graph using ggplot2
  p <- ggplot(data = plot.data, aes(x = X, y = rep(1:length(X)), color = Clusters)) +
    stat_ellipse() +
    ylab("Indexes of no relevance") +
    xlab("LDA1") +
    geom_point(aes(shape = Depth, size = 1.5)) +
    theme_bw()
  return(p)
}

# Import datasets
setwd("~/git/coralscape/results/")
struc.complex <- read.csv("~/git/coralscape/results/struc_complex_results.txt", sep = "\t")
#struc.complex <- read.csv("~/git/coralscape/results/axis_angle_rotation/struc_complex_results_X.csv")
#clusters <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/all-aga_1d_nc_20_6.csv")
clusters <- read.csv("all-aga_1d_wc_20_4_assigned.csv")
coordinates <- read.csv("~/git/coralscape/results/annotations.csv")

# Combine
metadata <- combine_metadata(struc.complex, clusters, coordinates)
rm(clusters, coordinates, struc.complex)

# Structure
str(metadata)
#metadata$Clusters[metadata$Clusters == 1] = "AA2"
#metadata$Clusters[metadata$Clusters == 2] = "AL2"
#metadata$Clusters[metadata$Clusters == 3] = "AH1"
#metadata$Clusters[metadata$Clusters == 4] = "AA1"
#metadata$Clusters[metadata$Clusters == 5] = "AL1"
#metadata$Clusters[metadata$Clusters == 6] = "AH2"
#metadata$Clusters[metadata$Clusters == 7] = "AH3"
              #AA1        #AA2        #AH1       #AH2      #AH3      #AL1      #AL2
colours <- c("#274e13","#8fce00", "#b06100", "#ffa500", "#ff4e00", "#6f0641", "#7223f0")
metadata[, c(2:4, 9)] <- data.frame(lapply(metadata[, c(2:4, 9)], as.factor))
str(metadata)
metadata$colony_rugosity_sqrt <- sqrt(metadata$colony_rugosity)
metadata$environment_rugosity_sqrt <- sqrt(metadata$environment_rugosity)
metadata$colony_elevation_corr <- fix_angles(metadata$colony_elevation)
metadata$ground_elevation_corr <- fix_angles(metadata$ground_elevation)

metadata$overhang_prop_2 <- NA
for (i in 1:length(metadata$overhang_prop)) {
  if (metadata$overhang_prop[i] > 1) {
    metadata$overhang_prop_2[i] <- 1
  } else {
    metadata$overhang_prop_2[i] <- metadata$overhang_prop[i]
  }
}
# AREA RESULTS DON'T MATCH 25/5/22 ####
#area_results <- read.delim("~/git/coralscape/results/area_results.txt")
#area_results_WP20 <- read.delim("~/git/coralscape/results/results_WP20/area_results_WP20.txt")
#area_results <- rbind(area_results, area_results_WP20)
#area_results_sub <-  subset(area_results, sample_name %in% metadata$Individual)
#area_results_sort <-  area_results_sub[order(area_results_sub[,2]),]
#identical(area_results_sort[,2], metadata[,1])
#metadata$area <- area_results_sort$colony_threeD_area
#metadata$area_sqrt <- sqrt(metadata$area)
#rm(area_results)
#rm(area_results_sort)
#rm(area_results_sub)

write_csv(metadata, file = "~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/metadata_X.csv")
# List of samples to remove due to issues ####
# See sample_investigation.xcls
env <- metadata[, c(18,24,26,27,28)]
env_x <- env
env_x$Clusters <- metadata$Clusters
env_x$Individual <- metadata$Individual
env_x$colony_points <- metadata$colony_points
to_remove <- NA
# Individuals that assigned to LM but are AC or HU
to_remove <- c("KP0529_LM_WP20", "KP0639_LM_SB20", "KP0399_LM_WP20", "KP0583_LM_WP20")
# Individuals that assigned to AC but are LM (or HU) 
to_remove <- append(to_remove, c("KP0765_AC_WP10", "KP0851_AC_WP10", "KP0108_AC_SB10", "KP0096_AC_SB10"))
# Indivudal that appears to be AA2, only next one along (771) but assigned to AA1, 590 AA2 but seems like AA1 need to check
to_remove <- append(to_remove, c("KP0769_AC_WP10", "KP0771_AC_WP10", "KP0590_AC_WP20"))
# Funky colony pics, the 650 lamarcki just some adjustment but mostly correct
to_remove <- append(to_remove, c( "KP0650_LM_SB20", "KP0774_AC_WP10", "KP0400_LM_WP20", "KP0350_LM_WP20"))
# Strangely clones with LM
to_remove <- append(to_remove, c("KP0604_AC_SB20"))
# Funky measures
#to_remove <- append(to_remove,(c("KP0293_LM_WP20",)))
metadata <- metadata[!metadata$Individual %in% to_remove,]
env_x <- env_x[!env_x$Individual %in% to_remove,]
env_x$area <- metadata$area
env_x$colony_rugosity_sqrt <- metadata$colony_rugosity_sqrt
env_x$colony_rugosity <- metadata$colony_rugosity
env_x$colony_elevation_corr <- metadata$colony_elevation_corr
env_x$overhang_prop_2 <- metadata$overhang_prop_2
env_x$environment_rugosity <- metadata$environment_rugosity
env_x$site <- metadata$Site

write.csv(env_x, "env_x.csv", quote = FALSE)
# Plots ####
k_boxplot(metadata, metadata$Clusters, metadata$environment_rugosity_sqrt, "Square-root of Environment Rugosity", colours)
k_boxplot(metadata, metadata$Clusters, metadata$environment_rugosity_sqrt, "Square-root of Environment Rugosity", colours, facet = "All")
k_boxplot(metadata, metadata$Clusters, metadata$environment_rugosity_sqrt, "Square-root of Envirosnment Rugosity", colours, facet = "Depth")
#k_boxplot(metadata, metadata$Clusters, metadata$colony_elevation, "Colony Elevation", colours)
#k_boxplot(metadata, metadata$Clusters, metadata$colony_elevation, "Colony Elevation", colours, facet = "All")
#k_boxplot(metadata, metadata$Clusters, metadata$colony_elevation, "Colony Elevation", colours, facet = "Depth")
k_boxplot(metadata, metadata$Clusters, metadata$colony_elevation_corr, "Colony Elevation corr", colours)
k_boxplot(metadata, metadata$Clusters, metadata$colony_elevation_corr, "Colony Elevation corr", colours, facet = "All")
k_boxplot(metadata, metadata$Clusters, metadata$colony_elevation_corr, "Colony Elevation corr", colours, facet = "Depth")

k_boxplot(metadata, metadata$Clusters, metadata$ground_elevation_corr, "Ground Elevation corr", colours)
k_boxplot(metadata, metadata$Clusters, metadata$ground_elevation_corr, "Ground Elevation corr", colours, facet = "All")
k_boxplot(metadata, metadata$Clusters, metadata$ground_elevation_corr, "Ground Elevation corr", colours, facet = "Depth")

k_boxplot(metadata, metadata$Clusters, metadata$colony_rugosity_sqrt, "Square-root of Colony Rugosity", colours)
k_boxplot(metadata, metadata$Clusters, metadata$colony_rugosity_sqrt, "Square-root of Colony Rugosity", colours, facet = "All")
k_boxplot(metadata, metadata$Clusters, metadata$colony_rugosity_sqrt, "Square-root of Colony Rugosity", colours, facet = "Depth")
k_boxplot(metadata, metadata$Clusters, metadata$colony_rugosity, "Colony Rugosity", colours, facet = "Depth")
#k_boxplot(metadata, metadata$Clusters, metadata$overhang_prop, "Overhang prop", colours)
k_boxplot(metadata, metadata$Clusters, metadata$overhang_prop_2, "Overhang proportion", colours)
k_boxplot(metadata, metadata$Clusters, metadata$overhang_prop_2, "Overhang proportion", colours, facet = "All")
k_boxplot(metadata, metadata$Clusters, metadata$overhang_prop_2, "Overhang proportion", colours, facet = "Depth")

k_boxplot(metadata, metadata$Clusters, metadata$outcrop_prop, "Outcrop proportion", colours)
k_boxplot(metadata, metadata$Clusters, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "All")
k_boxplot(metadata, metadata$Clusters, metadata$outcrop_prop, "Outcrop proportion", colours, facet = "Depth")


ggplot(metadata, aes(colony_rugosity, colony_elevation_corr, colour = Clusters)) +
  geom_point() + theme_classic()
  
param_list = c('colony_rugosity_sqrt', 'environment_rugosity_sqrt', 'colony_elevation_corr', 'outcrop_prop', 'overhang_prop_2', 'area_sqrt')
site_list = c('SB05', 'SB10', 'SB20', 'WP05', 'WP10', 'WP20', "CA05")
depth_list = c('20', '10', '5')

for (i in param_list) {
  print("Starting new parameter loop")
  p <- k_boxplot(metadata, metadata$Clusters, metadata[,i], i, colours)
  print(p)
  ggsave(paste0("plots/", i, ".pdf"), p)
  p1 <- k_boxplot(metadata, metadata$Clusters, metadata[,i], i, colours, facet = "All")
  print(p1)
  ggsave(paste0("plots/", i, "_facet.pdf"), p1, height = 30, width = 20, units = "cm")
  print("Starting site loop")
  pd <- k_boxplot(metadata, metadata$Clusters, metadata[,i], i, colours, facet = "Depth")
  print(pd)
  ggsave(paste0("plots/", i, "_depth.pdf"), pd, height = 15, width = 20, units = "cm")
  pl <- k_boxplot(metadata, metadata$Clusters, metadata[,i], i, colours, facet = "Loc")
  print(pl)
  ggsave(paste0("plots/", i, "_loc.pdf"), pl, height = 15, width = 20, units = "cm")
  for (j in site_list) {
    site <- metadata[metadata$Site == j,]
    p2 <- k_boxplot(site, site$Clusters, site[,i], paste(i, j), colours)
    print(p2)
    ggsave(paste0("plots/", i, "_", j, ".pdf"), p2)}
  print("Starting depth loop")
  for (k in depth_list) {
      depth <- metadata[metadata$Depth == k,]
      p3 <- k_boxplot(depth, depth$Clusters, depth[,i], paste(i, k), colours)
      print(p3)
      ggsave(paste0("plots/", i, "_", k, ".pdf"), p3)}
  }

# Proportion of colonies ###
table(metadata$Clusters, metadata$Depth)
ggplot(metadata, aes(Depth, fill = Clusters)) +
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

ggplot(metadata, aes(Clusters, fill = Depth)) +
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

# Other data format for elevation ####
metadata$elevation_group <- NA
#for (i in 1:nrow(metadata)) {
#  if (metadata$overhang_prop_2[i] >= 0.80 | metadata$colony_elevation[i] >= 90) {
#    metadata$elevation_group[i] = 7
#  } else if (metadata$overhang_prop_2[i] >= 0.60) {
#    metadata$elevation_group[i] = 6
#  } else if (metadata$overhang_prop_2[i] >= 0.40) {
#    metadata$elevation_group[i] = 5
#  } else if (metadata$overhang_prop_2[i] >= 0.20) {
#    metadata$elevation_group[i] = 4
#  } else if (metadata$colony_elevation[i] <= 30) {
#    metadata$elevation_group[i] = 1
#  } else if (metadata$colony_elevation[i] <= 60) {
#    metadata$elevation_group[i] = 2
#  } else if (metadata$colony_elevation[i] < 90) {
#    metadata$elevation_group[i] = 3
#  } else {
#    print('something wrong')
#  }
#}
for (i in 1:nrow(metadata)) {
if (metadata$overhang_prop_2[i] >= 0.80) {
  metadata$elevation_group[i] = 7
} else if (metadata$overhang_prop_2[i] >= 0.60) {
  metadata$elevation_group[i] = 6
} else if (metadata$overhang_prop_2[i] >= 0.40) {
  metadata$elevation_group[i] = 5
} else if (metadata$overhang_prop_2[i] >= 0.20) {
  metadata$elevation_group[i] = 4
} else if (metadata$colony_elevation_corr[i] <= 30) {
  metadata$elevation_group[i] = 1
} else if (metadata$colony_elevation_corr[i] <= 60) {
  metadata$elevation_group[i] = 2
} else if (metadata$colony_elevation_corr[i] <= 90) {
  metadata$elevation_group[i] = 3
} else {
  print('something wrong')
}}


metadata$elevation_group
boxplot(metadata$elevation_group~metadata$Clusters)
k_boxplot(metadata, metadata$Clusters, metadata$elevation_group, "Elevation group", colours)
k_boxplot(metadata, metadata$Clusters, metadata$elevation_group, "Elevation group", colours, facet = "All")

metadata$elevation_group <- as.factor(metadata$elevation_group)
ggplot(metadata, aes(Clusters, fill = elevation_group)) +
  geom_bar(position = "fill") +
  scale_y_continuous(label = scales::percent) +
  scale_fill_manual(values = c(blues9[c(3,5,7)], "lightgrey", "grey", "darkgrey", "black")) +
  theme(axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18, face = 'bold')) +
  ylab("Proportion") +
  ggtitle("Cluster vs. Elevation group") +
  theme_classic()

ggplot(metadata, aes(Clusters, fill = elevation_group)) +
  geom_bar(position = "fill") +
  scale_y_continuous(label = scales::percent) +
  scale_fill_manual(values = c(blues9[c(3,5,7)], "lightgrey", "grey", "darkgrey", "black")) +
  theme(axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18, face = 'bold')) +
  ylab("Proportion") +
  ggtitle("Cluster vs. Elevation group") +
  theme_classic() +
  facet_wrap(~ Depth) +
  coord_flip()

library(gmodels)
CrossTable(metadata$Clusters, metadata$elevation_group,
           format = 'SPSS', expected = T, prop.chisq = T)

# Does size really matter? ####
boxplot(metadata$area ~ metadata$Clusters)
k_boxplot(metadata, metadata$Clusters, metadata$area, "3D colony surface area", colours)
k_boxplot(metadata, metadata$Clusters, metadata$area, "3D colony surface area", colours, facet = "Yes")
k_boxplot(metadata, metadata$Clusters, metadata$area_sqrt, "Square root of 3D colony surface area", colours)
a <- k_boxplot(metadata, metadata$Clusters, metadata$area_sqrt, "Square root of 3D colony surface area", colours)
a
ggsave(paste0("plots/area_clust.png"), a)
a1 <- k_boxplot(metadata, metadata$Clusters, metadata$area_sqrt, "Square root of 3D colony surface area ", colours, facet = "All")
a1
ggsave(paste0("plots/area_clust_facet.png"), a1, height = 30, width = 20, units = "cm")
a2 <- k_boxplot(metadata, metadata$Clusters, metadata$area_sqrt, "Square root of 3D colony surface area ", colours, facet = "Depth")
a2
ggsave(paste0("plots/area_clust_depth.png"), a2, height = 30, width = 20, units = "cm")
# Does size really matter? - Rugosity ####
ggplot(metadata, aes(area_sqrt, colony_rugosity_sqrt, colour = Clusters)) +
  geom_point() +
  scale_color_manual(name = "Clusters", values = colours) +
  theme_classic(base_size = 12)

ggplot(metadata, aes(area_sqrt, colony_rugosity_sqrt, colour = Clusters)) +
  geom_point() +
  scale_color_manual(name = "Clusters", values = colours) +
  theme_classic(base_size = 12) +
  geom_smooth(method = "lm", se = T, alpha = 0.15)

ggplot(metadata, aes(area_sqrt, colony_rugosity_sqrt, colour = Clusters)) +
  geom_point() +
  scale_color_manual(name = "Clusters", values = colours) +
  theme_classic(base_size = 12) +
  facet_wrap(~ Depth + Loc)

ggplot(metadata, aes(area_sqrt, colony_rugosity_sqrt, colour = Clusters)) +
  geom_point() +
  scale_color_manual(name = "Clusters", values = colours) +
  theme_classic(base_size = 12) +
  facet_wrap(~ Depth + Loc) +
  geom_smooth(method = "lm", se = T, alpha = 0.15)


# Does size really matter? - Colony elevation ####

ggplot(metadata, aes(area_sqrt, colony_elevation, colour = Clusters)) +
  geom_point() +
  scale_color_manual(name = "Clusters", values = colours) +
  theme_classic(base_size = 12)

ggplot(metadata, aes(area_sqrt, colony_elevation, colour = Clusters)) +
  geom_point() +
  scale_color_manual(name = "Clusters", values = colours) +
  theme_classic(base_size = 12) +
  geom_smooth(method = "lm", se = T, alpha = 0.15)

ggplot(metadata, aes(area_sqrt, colony_elevation, colour = Clusters)) +
  geom_point() +
  scale_color_manual(name = "Clusters", values = colours) +
  theme_classic(base_size = 12) +
  facet_wrap(~ Depth + Loc)

ggplot(metadata, aes(area_sqrt, colony_elevation, colour = Clusters)) +
  geom_point() +
  scale_color_manual(name = "Clusters", values = colours) +
  theme_classic(base_size = 12) +
  facet_wrap(~ Depth + Loc) +
  geom_smooth(method = "lm", se = T, alpha = 0.15)

# *** Statistics *** ####

# Checking normality
shapiro.test(metadata$colony_elevation) # not-normal
shapiro.test(metadata$colony_rugosity_sqrt) # not normal, from 1 to inf
shapiro.test(metadata$environment_rugosity_sqrt) #normal
shapiro.test(metadata$outcrop_prop) # not normal
shapiro.test(metadata$overhang_prop) # not normal, from 0 to 1, zero-inflated
#for (i in (levels(metadata$Clusters))) {
#  hist(metadata[metadata$Clusters == i,]$colony_elevation, main = paste(i, "elevation"))
#  hist(metadata[metadata$Clusters == i,]$colony_rugosity, main = paste(i, "rugosity"))
#  hist(metadata[metadata$Clusters == i,]$overhang_prop, main = paste(i, "overhang"))
#  hist(metadata[metadata$Clusters == i,]$outcrop_prop, main = paste(i, "outcrop prop"))
#  hist(metadata[metadata$Clusters == i,]$environment_rugosity, main = paste(i, "env rugosity"))
#}

# Elevation ####
lm.colony_elevation <- lm(colony_elevation_corr ~ Clusters, data = metadata)
summary(lm.colony_elevation)
TukeyHSD(aov(lm.colony_elevation))

lmer.colony_elevation <- lmer(colony_elevation_corr ~ Clusters + (1|Loc), data = metadata)
summary(lmer.colony_elevation)
emmeans(lmer.colony_elevation, list(pairwise ~ Clusters), adjust = "tukey")

lmer.colony_elevation.depth <- lmer(colony_elevation_corr ~ Depth + (1|Loc), data = metadata)
summary(lmer.colony_elevation.depth)
emmeans(lmer.colony_elevation.depth, list(pairwise ~ Depth), adjust = "tukey") # No differences between depths?

lmer2.colony_elevation <- lmer(colony_elevation_corr ~ Clusters*Depth + (1|Loc), data = metadata)
summary(lmer2.colony_elevation)
emmeans(lmer2.colony_elevation, list(pairwise ~ Clusters:Depth), adjust = "tukey") # No interactions

lmer3.colony_elevation <- lmer(colony_elevation_corr ~ Clusters + (1|Loc/Depth), data = metadata)
summary(lmer3.colony_elevation)
emmeans(lmer3.colony_elevation, list(pairwise ~ Clusters), adjust = "tukey")

lmer4.colony_elevation <- lmer(colony_elevation_corr ~ Clusters + (Clusters|Loc/Depth), data = metadata)
summary(lmer4.colony_elevation)
# NaNs issues...

# tests with difference variance...
kruskal.test(metadata$colony_elevation_corr, metadata$Clusters)
boxplot(metadata$colony_elevation_corr ~ metadata$Clusters)
library(FSA)
dunnTest(metadata$colony_elevation ~ metadata$Clusters, 
         method = 'bonferroni')

# Elevation - Testing with just AA ####
results.colony_elevation.aa <- data.frame()
for (i in 1:100) {
metadata.aa1 <- metadata[metadata$Clusters == "AA1",]
metadata.aa2 <- metadata[metadata$Clusters == "AA2",]
metadata.aa1.sample <- metadata.aa1[sample(1:nrow(metadata.aa1), length(metadata.aa2$Clusters)),]
metadata.aa <- rbind(metadata.aa1.sample, metadata.aa2)

lm.colony_elevation.aa <- lm(colony_elevation ~ Clusters, data = metadata.aa)
summary(lm.colony_elevation.aa)
TukeyHSD(aov(lm.colony_elevation.aa))$Clusters
row <- TukeyHSD(aov(lm.colony_elevation.aa))$Clusters
results.colony_elevation.aa <- rbind(results.colony_elevation.aa, row)

lmer.colony_elevation.aa <- lmer(colony_elevation ~ Clusters + (1|Loc), data = metadata.aa)
summary(lmer.colony_elevation.aa)
emmeans(lmer.colony_elevation.aa, list(pairwise ~ Clusters), adjust = "tukey")

row2 <- emmeans(lmer.colony_elevation.aa, list(pairwise ~ Clusters), adjust = "tukey")$`pairwise differences of Clusters`
print(row2)
} # 60% signifcant


# Environment rugoisty ####
lm.environment_rugosity <- lm(environment_rugosity_sqrt ~ Clusters, data = metadata)
summary(lm.environment_rugosity)
TukeyHSD(aov(lm.environment_rugosity))
plot(lm.environment_rugosity)

lmer.environment_rugosity <- lmer(environment_rugosity_sqrt ~ Clusters + (1|Loc/Depth), data = metadata)
summary(lmer.environment_rugosity)
emmeans(lmer.environment_rugosity, list(pairwise ~ Clusters), adjust = "tukey")

lmer.environment_rugosity.depth <- lmer(environment_rugosity_sqrt ~ Depth + (1|Loc), data = metadata)
summary(lmer.environment_rugosity.depth)
emmeans(lmer.environment_rugosity.depth, list(pairwise ~ Depth), adjust = "tukey") # SUPER DIFF WITH DEPTH

lmer2.environment_rugosity <- lmer(environment_rugosity_sqrt ~ Clusters * Depth + (1|Loc), data = metadata)
summary(lmer2.environment_rugosity)
emmeans(lmer2.environment_rugosity, list(pairwise ~ Clusters:Depth), adjust = "tukey")
# Doesn't neccessary have to do with the clusters just that environment rugosity changes with depth

# Outcrop prop ####
lm.outcrop_prop <- lm(outcrop_prop ~ Clusters, data = metadata)
summary(lm.outcrop_prop)
TukeyHSD(aov(lm.outcrop_prop))

lm.outcrop.depth <- lm(outcrop_prop ~ Depth, data = metadata)
summary(lm.outcrop.depth)
TukeyHSD(aov(lm.outcrop.depth))

lmer.outcrop_prop <- lmer(outcrop_prop ~ Clusters*Depth + (1|Loc), data = metadata)
summary(lmer.outcrop_prop)
emmeans(lmer.outcrop_prop, list(pairwise ~ Clusters:Depth), adjust = "tukey")
#AA2 Depth10 - AA2 Depth20 -0.05981 0.0161 348.6  -3.705  0.0352

lm.outcrop_prop.depth <- lm(outcrop_prop ~ Clusters*Depth, data = metadata)
summary(lm.outcrop_prop.depth)
emmeans(lm.outcrop_prop.depth, list(pairwise ~Clusters*Depth), adjust = "tukey")

glmer.outcrop.interaction <- glmer(outcrop_prop ~ Clusters*Depth + (1|Loc), family = 'binomial', data = metadata)
emmeans(glmer.outcrop.interaction, list(pairwise~Clusters*Depth))

# Incorporating an interactivate effect with depth
lm.outcrop.interaction <- lm(outcrop_prop ~ Clusters:Depth, data = metadata)
summary(lm.outcrop.interaction) 
anova(lm.outcrop.interaction)
TukeyHSD(aov(lm.outcrop.interaction))
plot(lm.outcrop.interaction)

lmer.outcrop_prop <- lmer(outcrop_prop ~ Clusters + (1|Loc/Depth), data = metadata)
summary(lmer.outcrop_prop) 
emmeans(lmer.outcrop_prop, list(pairwise~Clusters), adjust = "tukey")
lmer.outcrop_prop2 <- lmer(outcrop_prop ~ Clusters + (Clusters|Loc/Depth), data = metadata)
summary(lmer.outcrop_prop2) 
emmeans(lmer.outcrop_prop2, list(pairwise~Clusters), adjust = "tukey")

# Overhang ####
lm.overhang <- lm(overhang_prop ~ Clusters, data = metadata)
summary(lm.overhang)
anova(lm.overhang)
TukeyHSD(aov(lm.overhang))

lmer.overhang <- lmer(overhang_prop ~ Clusters + (1|Loc), data = metadata)
summary(lmer.overhang)
emmeans(lmer.overhang, list(pairwise~Clusters), adjust = "tukey")

lmer2.overhang <- lmer(overhang_prop ~ Clusters + (1|Loc/Depth), data = metadata)
summary(lmer2.overhang)
emmeans(lmer2.overhang, list(pairwise~Clusters), adjust = "tukey")


# Colony rugosity ####
lm.colony_rugosity <- lm(colony_rugosity_sqrt ~ Clusters, data = metadata)
summary(lm.colony_rugosity)
TukeyHSD(aov(lm.colony_rugosity))
plot(lm.colony_rugosity)

lm.colony_rugosity.20 <- lm(colony_rugosity_sqrt ~ Clusters, data = metadata[metadata$Depth == "20",])
summary(lm.colony_rugosity)
TukeyHSD(aov(lm.colony_rugosity))

lmer.colony_rugosity <- lmer(colony_rugosity_sqrt ~ Clusters + (1|Loc/Depth), data = metadata)
emmeans(lmer.colony_rugosity, list(pairwise ~ Clusters), adjust = "tukey")

lmer.colony_rugosity <- lmer(colony_rugosity_sqrt ~ Clusters + (1|Loc), data = metadata)
emmeans(lmer.colony_rugosity, list(pairwise ~ Clusters), adjust = "tukey")

# Colony rugosity - Testing with just AA ####
results.colony_rugosity.aa <- data.frame()
for (i in 1:100) {
  metadata.aa1 <- metadata[metadata$Clusters == "AA1",]
  metadata.aa2 <- metadata[metadata$Clusters == "AA2",]
  metadata.aa1.sample <- metadata.aa1[sample(1:nrow(metadata.aa1), length(metadata.aa2$Clusters)),]
  metadata.aa <- rbind(metadata.aa1.sample, metadata.aa2)
  
  lm.colony_rugosity.aa <- lm(colony_rugosity ~ Clusters, data = metadata.aa)
  summary(lm.colony_rugosity.aa)
  TukeyHSD(aov(lm.colony_rugosity.aa))$Clusters
  row <- TukeyHSD(aov(lm.colony_rugosity.aa))$Clusters
  results.colony_rugosity.aa <- rbind(results.colony_rugosity.aa, row)
  
  lmer.colony_rugosity.aa <- lmer(colony_rugosity ~ Clusters + (1|Loc), data = metadata.aa)
  summary(lmer.colony_rugosity.aa)
  emmeans(lmer.colony_rugosity.aa, list(pairwise ~ Clusters), adjust = "tukey")
  
  row2 <- emmeans(lmer.colony_rugosity.aa, list(pairwise ~ Clusters), adjust = "tukey")$`pairwise differences of Clusters`
  print(row2)
} # not significant only 3/100


# Mixed-models ####
mm.elevation <- lm(colony_elevation_corr ~ Clusters*colony_rugosity, data = metadata)
summary(mm.elevation)
TukeyHSD(aov(mm.elevation))

mmer.elevation <- lmer(colony_elevation_corr ~ Clusters*colony_rugosity + (1|Loc/Depth), data = metadata)
summary(mmer.elevation)
emmeans(mmer.elevation, list(pairwise ~ Clusters:colony_rugosity), adjust = "tukey")

# Scaling ####
env <- metadata[, c(24,32,33,34)]
env.scale <- env
env.scale$colony_elevation_corr <- scale(env$colony_elevation_corr)
env.scale$overhang_prop_2 <- scale(env$overhang_prop_2)
env.scale$outcrop_prop <- scale(env$outcrop_prop)
env.scale$environment_rugosity_sqrt <- scale(env$environment_rugosity_sqrt)
env.scale$z <- scale(metadata$z)
#env.scale$area_sqrt <- scale(env$area_sqrt)
env.scale$Clusters <- metadata$Clusters

# PCA ####
# Visualise data
pairs(env.scale, lower.panel = NULL, col = as.numeric(env.scale$Clusters))

pca <- prcomp(env.scale[,-5], scale = T)
biplot(pca)

# PCA with rda function ####
rda <- rda(env.scale[,-5], scale = T)
biplot(rda)

#make basic plot
cluster.names <- levels(env.scale$Clusters)
biplot(rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))

#Add "hulls"
ordihull(rda,
         group = env.scale$Clusters, col = colours)
legend("topleft",
       col = colours,
       lty = 1,
       legend = cluster.names)
title("PCA of coral taxa based on microhabitat variables")
with(env.scale, levels(Clusters))
colvec <- colours
with(env.scale, points(rda, display = "sites", col = colvec[Clusters], pch = 21, bg = colvec[Clusters]))

#pca.centroids <- aggregate(rda$CA$u, list(Type = env.scale$Clusters), mean)
#points(pca.centroids$PC1, pca.centroids$PC2, col = colvec)

dis <- vegdist(rda$CA$Xbar, "euclidean")
centroid <- betadisper(dis, env.scale$Clusters, "centroid")

# Different depths #####
k_rda(env.scale, "5")
k_rda(env.scale, "10")
k_rda(env.scale, "20")
# TODO: colours for RDA


# Linear Discriminant Analysis ####
# LDA is like PCA but it focuses on maximising the seperability
# among known categories
env.scale.lda <- env.scale[,-5]
env.scale.lda$Clusters <- metadata$Clusters
lda <- lda(formula = Clusters ~ ., data = env.scale)

## get the x,y coordinates for the LDA plot
data.lda.values <- predict(lda)
env.scale.lda$Depth <- metadata$Depth
#env.scale.lda <- env.scale.lda %>% unite(cluster.depth, Clusters, Depth)
env.scale.lda$Clusters <- metadata$Clusters
## create a dataframe that has all the info we need to draw a graph
plot.data <- data.frame(LD1 = data.lda.values$x[,1], LD2 = data.lda.values$x[,2], LD3 =  data.lda.values$x[,3],
                        LD4 = data.lda.values$x[,4], LD5 = data.lda.values$x[,5],
                        Clusters = env.scale.lda$Clusters, Depth = env.scale.lda$Depth)

head(plot.data)

## draw a graph using ggplot2
p <- ggplot(data = plot.data, aes(x = LD1, y = LD2, color = Clusters)) +
  scale_color_manual(values = colours) +
  theme_bw() +
  stat_ellipse() +
  ggtitle("LDA") +
  geom_point(aes(shape = Depth, alpha = 0.2))
p

p <- ggplot(data = plot.data, aes(x = LD3, y = LD4, color = Clusters)) +
  scale_color_manual(values = colours) +
  theme_bw() +
  stat_ellipse() +
  ggtitle("LDA") +
  geom_point(aes(shape = Depth, alpha = 0.2))
p


# LDA - for clusters ####
k_lda(env.scale.lda, "AA")
k_lda(env.scale.lda, "AH")
k_lda(env.scale.lda, "AL")

# LDA - Subsmaple AA1 of 20 m then vs. AA2 ####
env.scale.aa1 <- env.scale.lda[env.scale.lda$Clusters == "AA1" & env.scale.lda$Depth == "20",]
env.scale.aa2 <- env.scale.lda[env.scale.lda$Clusters == "AA2" & env.scale.lda$Depth == "20",]
env.scale.aa1.sample <- env.scale.aa1[sample(1:nrow(env.scale.aa1), length(env.scale.aa2$Clusters)),]
env.scale.aa <- rbind(env.scale.aa1.sample, env.scale.aa2)
env.scale.aa.lda <- env.scale.aa[,-7]
env.scale.aa.lda$Clusters <- droplevels(env.scale.aa.lda$Clusters)
lda <- lda(formula = Clusters ~ ., data = env.scale.aa.lda)

## get the x,y coordinates for the LDA plot
data.lda.values <- predict(lda)

## create a dataframe that has all the info we need to draw a graph
plot.data <- data.frame(X = data.lda.values$x[,1], Clusters = env.scale.aa$Clusters, Depth = env.scale.aa$Depth)

head(plot.data)

## draw a graph using ggplot2
p <- ggplot(data = plot.data, aes(x = X, y = rep(1:length(X)), color = Clusters)) +
  stat_ellipse() +
  xlab("LDA1") +
  ylab("Indexes of no relevance") +
  geom_point(aes(shape = Depth, size = 1.5)) +
  theme_bw()
p


# NMDS ####
rugosity <- aggregate(env.scale$environment_rugosity, list(Clusters = env.scale$Clusters), mean)
colnames(rugosity) <- c("Clusters", "env_rugosity")
elevation <- aggregate(env.scale$colony_elevation, list(Clusters = env.scale$Clusters), mean)
colnames(elevation) <- c("Clusters", "elevation")
outcrop <- aggregate(env.scale$outcrop_prop, list(Clusters = env.scale$Clusters), mean)
colnames(outcrop) <- c("Clusters", "outcrop")
overhang <- aggregate(env.scale$overhang_prop, list(Clusters = env.scale$Clusters), mean)
colnames(overhang) <- c("Clusters", "overhang")

mean_taxa <- cbind(rugosity, elevation$elevation, outcrop$outcrop, overhang$overhang)
colnames(mean_taxa) <- c("Clusters", "rugosity", "elevation", "outcrop", "overhang")
row.names(mean_taxa) <- mean_taxa[,1]
mean_taxa <- mean_taxa[,-1]
head(mean_taxa)
cluster.dist <- dist(mean_taxa, "euclidean")

#nmds <- metaMDSiter(cluster.dist)
nmds <- metaMDS(cluster.dist)
ordiplot(nmds, type = 'text')
envfit <- envfit(nmds, mean_taxa)
plot(envfit)

# Within species variation ####
species = "AA1"
ggplot(metadata[metadata$Clusters == species,], aes(overhang_prop_2, colony_elevation_corr, colour = Depth)) +
  geom_point() +
  scale_color_manual(name = "Depth", values = blues9[c(3, 6, 9)]) +
  theme_classic(base_size = 12)
 
ggplot(metadata[metadata$Clusters == species,], aes(outcrop_prop, colony_elevation_corr, colour = Depth)) +
  geom_point() +
  scale_color_manual(name = "Depth", values = blues9[c(3, 6, 9)]) +
  theme_classic(base_size = 12)

ggplot(metadata[metadata$Clusters == species,], aes(environment_rugosity_sqrt, colony_elevation_corr, colour = Depth)) +
  geom_point() +
  scale_color_manual(name = "Depth", values = blues9[c(3, 6, 9)]) +
  theme_classic(base_size = 12)

ggplot(metadata[metadata$Clusters == species,], aes(area_sqrt, colony_elevation_corr, colour = Depth)) +
  geom_point() +
  scale_color_manual(name = "Depth", values = blues9[c(3, 6, 9)]) +
  theme_classic(base_size = 12)
# if need  +
#geom_smooth(method = "lm", se = T, alpha = 0.15)

ggplot(metadata[metadata$Clusters == species,], aes(colony_elevation_corr, Depth)) +
  geom_boxplot()
ggplot(metadata[metadata$Clusters == species,], aes(colony_rugosity_sqrt, colony_elevation_corr, color = Depth)) +
  geom_point() +
  theme_classic(base_size = 12)

# AC growth forms####
subset_ac <- metadata[metadata$Clusters == "AA1" | metadata$Clusters == "AA2",]
subset_ac_20 <- subset_ac[subset_ac$Site == "SB20" | subset_ac$Site == "WP20",]
ac_growth_forms <- read.csv("~/git/coralscape/ac_growth_forms.csv", header = FALSE)
identical(ac_growth_forms[,1], subset_ac_20[,1])
subset_ac_20$growth <- ac_growth_forms$V2
subset_ac_20 <- subset_ac_20 %>% unite(growth_cluster, Clusters, growth, remove = FALSE)

ggplot(subset_ac_20, aes(z, growth, colour = Clusters)) +
  geom_boxplot() + theme_classic() + facet_wrap(~ Site)

ggplot(subset_ac_20, aes(colony_rugosity_sqrt, colony_elevation_corr, colour = growth_cluster)) +
  geom_boxplot() +
  theme_classic(base_size = 12)

ggplot(subset_ac_20, aes(colony_rugosity_sqrt, growth, colour = Clusters)) +
  geom_boxplot() +
  geom_point() +
  scale_color_manual(values = c(""))
  ylab("Colony growth from") +
  xlab("Square root of colony rugosity") +
  theme_classic(base_size = 12)

ggplot(subset_ac_20, aes(colony_elevation_corr, colony_rugosity_sqrt, colour = growth)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_classic()

ggplot(subset_ac_20, aes(colony_elevation_corr, z, colour = Clusters)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_classic()

ggplot(subset_ac_20, aes(outcrop_prop, z, colour = Clusters)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_classic()

ggplot(subset_ac_20, aes(outcrop_prop, growth_cluster, colour = Clusters)) +
  geom_boxplot() + theme_classic()

ggplot(subset_ac_20, aes(colony_rugosity_sqrt, z, colour = Clusters)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_classic()

ggplot(subset_ac_20, aes(overhang_prop_2, z, colour = Clusters)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_classic()

lm_elevation_ac <- lm(colony_elevation_corr ~ growth_cluster, data = subset_ac_20)
summary(lm_elevation_ac)
TukeyHSD(aov(lm_elevation_ac))
#AA1_bi/uni-AA1_bi     -22.727444 -40.47946 -4.975426 0.0043327
#AA1_uni-AA1_bi        -25.320715 -43.07273 -7.568697 0.0009857
lmer_elevation_ac <- lmer(colony_elevation_corr ~ growth_cluster + (1|Loc), data = subset_ac_20)
emmeans(lmer_elevation_ac, list(pairwise ~ growth_cluster), adjust = "tukey")

##### lm subset
subset_lm <- metadata[metadata$Clusters == "AL1" | metadata$Clusters == "AL2",]
subset_lm_20 <- subset_lm[subset_lm$Site == "SB20" | subset_lm$Site == "WP20",]
ggplot(subset_lm_20, aes(overhang_prop_2, colony_elevation_corr, colour = Clusters)) +
  geom_point() + theme_classic()

ggplot(subset_lm_20, aes(overhang_prop_2, area_sqrt, colour = Clusters)) +
  geom_point() + theme_classic()

ggplot(subset_lm_20, aes(colony_elevation_corr, area_sqrt, colour = Clusters)) +
  geom_point() + theme_classic()
