# Title: The effect of microniches on snp data using ordinations
# Author: Katharine Prata
# Date created: 10/01/21

# Packages
library(vcfR)
library(adegenet)
library(ggplot2)
library(vegan)
library(adespatial)

# Functions
combine_metadata <- function(genind, struc.complex, clusters, coordinates) {
  # match
  struc.complex.subset <- subset(struc.complex, X %in% indNames(genind))
  clusters.subset <- subset(clusters, Individual %in% indNames(genind))
  coordinates.subset <- subset(coordinates, X %in% indNames(genind))
  # sort
  sorted.clusters.subset <- clusters.subset[order(clusters.subset$Individual),]
  # check if identiical
  cat("Checking gen imp and clusters: ", identical(indNames(genind),sorted.clusters.subset[,1]))
  cat("\nChecking gen imp and struc complex: ", identical(indNames(genind),struc.complex.subset[,1]))
  cat("\nChecking gen imp and coordinates: ", identical(indNames(genind), coordinates.subset[,1]))
  # Put metadata together
  metadata <- sorted.clusters.subset
  metadata[,20:24] <- struc.complex.subset[,2:6]
  metadata[, 25:27] <- coordinates.subset[,2:4]
  return(metadata)
}

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

overhang_barplot <- function(metadata, factor) {
  overhang <- table(metadata$overhang, factor) %>% matrix(nrow = 2, ncol = length(levels(factor)))
  colnames(overhang) <- levels(factor)
  rownames(overhang) <- levels(metadata$overhang)
  barplot(overhang, main = "Presence of overhang", xlab = "Taxa",
          col = c("Red", "Green"), ylim = c(0,100))
  legend("topright", c("No", "Yes"), fill = c("Red", "Green"))
}

boxplot <- function(df, cat, cont) {
  
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
struc.complex <- read.csv("~/git/coralscape_open3d/results/sample_metadata.csv")
aga.genind <- vcfR2genind(read.vcfR("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/all-aga_1d_nc_20.vcf"))
aga.clusters <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/all-aga_1d_nc_20_6.csv")
aga.coordinates <- read.csv("~/git/coralscape_open3d/results/transformed_annotations.csv")

# Organise data
gen.subset <- subset(aga.genind, indNames(aga.genind) %in% struc.complex$X)
metadata <- combine_metadata(gen.subset, struc.complex, aga.clusters, aga.coordinates)
gen.imp <- organise_snp_matrix(gen.subset) # WARNING IMPUTATION METHOD PROBABLY INCORRECT NEED POPs

# genetic distance data
gen.subset$pop <- as.factor(indNames(gen.subset))
gen.pop <- genind2genpop(gen.subset)
distgenEUCL <- dist(gen.pop, method = "euclidean", diag = FALSE, upper = FALSE)
# The advantage of PCoA is that it doesn't take into account the missing observations
Pcoa <- ape::pcoa(distgenEUCL)
Pcoa
X <- Pcoa$vectors

# Clean working space
rm(aga.genind)
rm(struc.complex)
rm(aga.clusters)
rm(aga.coordinates)
rm(gen.pop)
rm(gen.subset)

# Data structure
str(metadata)
metadata[, c(2:6, 13, 24)] <- data.frame(lapply(metadata[, c(2:6, 13, 24)], as.factor))
levels(metadata$Clusters) <- c("AA1", "AL1", "AA2", "AL2") # have to change if using different datasets
str(metadata)

# choosing steepest angle
metadata$theta <- NA
for (row in 1:length(metadata$xz)) {
  if (metadata$xz[row] > metadata$yz[row]) {
    metadata$theta[row] <- metadata$xz[row]
  } else {
    metadata$theta[row] <- metadata$yz[row]
  }
}

# Environmental data ####
# Overhang
overhang_barplot(metadata, metadata$Species)  # Species
overhang_barplot(metadata, metadata$Clusters)  # Cluster

# Theta
boxplot(metadata, metadata$Species, metadata$theta) + ggtitle('Theta')
boxplot(metadata, metadata$Clusters, metadata$theta) + ggtitle('Theta')
# maybe colony rugosity would be good to include here!

# Relative depth proportion species
boxplot(metadata, metadata$Species, metadata$prop) + ggtitle('Outcrop position')
boxplot(metadata, metadata$Clusters, metadata$prop) + ggtitle('Outcrop position')

# Z
boxplot(metadata, metadata$Species, metadata$z) + ggtitle('Relative depth')
boxplot(metadata, metadata$Clusters, metadata$z) + ggtitle('Relative depth')

#### Linear Models ####
# glmer for hierarchy 
metadata$overhang.bin <- ifelse(metadata$overhang == "Yes", 1, 0)
glm.overhang <- glm(overhang ~ Clusters, data = metadata, family = binomial) # (1|Loc/Depth)
summary(glm.overhang)
# probability space with confidence intervals
emmeans()

lm.theta <- lm(theta ~ Clusters, data = metadata)
summary(lm.theta)
TukeyHSD(lm.theta)

#### RDA ####
full.rda <- rda(gen.imp ~ theta * overhang * prop * z, data = metadata, scale = T)
RsquareAdj(full.rda) # 0.09%
summary(eigenvals(full.rda, model = "constrained"))
#RDA1 explains 62% of the variations, RDA2 explains 6%. RDA3 4.4%, RDA4 3.9%
summary(eigenvals(full.rda, model = "unconstrained"))
# PC1 explains 43% of the variation, PC2 4%, PC3 0.3%, PC4 2%
summary(full.rda)

#Partitioning of correlations:
#Inertia Proportion
#Total           21764     1.0000
#Constrained      6062     0.2786
#Unconstrained   15702     0.7214

step(full.rda)
# removing each facotr the AIC decreases and removing clusters the AIC increases
# best model: Call: rda(formula = gen.subset.imp ~ overhang + z, data = metadata, scale = T)

# Significance
#signif.full.c <- anova.cca(full.rda)
#signif.full.c 
# Model is signicant 0.004

# Reduced rda ####
reduced.rda <- rda(gen.imp ~ overhang + z, data = metadata, scale = T)
RsquareAdj(reduced.rda) # 9%
summary(eigenvals(reduced.rda, model = "constrained"))
summary(eigenvals(reduced.rda, model = "unconstrained"))
summary(reduced.rda)

#Partitioning of correlations:
#Inertia Proportion
#Total           21764     1.0000
#Constrained      2696     0.1239
#Unconstrained   19068     0.8761

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
#  > how <- how(nperm=499, plots = Plots(strata=metadata$Loc))
# > anova(ord, by="term", permutations = how)

# Do a LDA




# Do a dbMEM ####
coords <- metadata[, 25:27]
dist.coords <- dist(coords)
dbMEM <- dbmem(dist.coords, MEM.autocor = c('positive')) # thresh chosen by the length of the lingest edge of the minimum spanning tree
dbMEM
attributes(dbMEM)$values # postive eigenvalues
attr(dbMEM, "listw")
adegraphics::s.label(coords[,1:2], attr(dbMEM, "listw")) # specify where to find function
s.value(coords[,1:2], dbMEM[,18])
plot(rownames(dbMEM), dbMEM$MEM1) # sine wave

Y <- cbind(dbMEM[,1:18], metadata$theta, metadata$overhang, metadata$prop, metadata$z)
cor(Y[,-20]) #  0.479627128 z and prop correlated and mems and other values
# what should be the minimum correlation

# dbRDA ####
dbMEM.rda <- dbrda(X ~ ., data = Y, scale = T)
vif.cca(dbMEM.rda) # When VIF is greater than around 10, this is problematic.

RsquareAdj(dbMEM.rda) # 13%
summary(eigenvals(dbMEM.rda, model = "constrained"))
summary(eigenvals(dbMEM.rda, model = "unconstrained"))
summary(dbMEM.rda)

#Partitioning of mean squared Euclidean distance:
#Inertia Proportion
#Total           13339     1.0000
#Constrained      5271     0.3952
#Unconstrained    8068     0.6048

# Significance
signif.dbMEM.rda <- anova.cca(dbMEM.rda)
signif.dbMEM.rda  # significant

signif.term.dbMEM.rda <- anova.cca(dbMEM.rda, by = 'terms')
signif.term.dbMEM.rda # MEM12 significant, overhang, z and nearly prop

dbMEM.rda0 <- capscale(X ~ 1, Y)
Sel <- ordiR2step(dbMEM.rda0, scope = formula(dbMEM.rda), direction = 'both')
Sel$anova #only overhang and prop
ade4::s.value(coords[,1:2], dbMEM[,12]) # d = 5?

# Select terms
Ysel <- as.data.frame(cbind(dbMEM$MEM12, metadata$overhang, metadata$z))
colnames(Ysel) <- c("MEM12", "Overhang", "prop")
dbMEM.rdaS <- dbrda(X ~ ., data = Ysel)
summary(dbMEM.rdaS)
RsquareAdj(dbMEM.rdaS) # 13%
signif.dbMEM.rdaS <- anova.cca(dbMEM.rdaS)
signif.dbMEM.rdaS # significant
signif.term.dbMEM.rdaS <- anova.cca(dbMEM.rdaS, by = 'terms')
signif.term.dbMEM.rdaS # all significant!

site <- cbind(scores(dbMEM.rdaS, display = "sites", choices = c(1,2,3,4), scaling = 1), metadata$Clusters)
AA1 <- site[site[,5] == 1,]
AL1 <- site[site[,5] == 2,]
AA2 <- site[site[,5] == 3,]
AL2 <- site[site[,5] == 4,]

summary(eigenvals(dbMEM.rdaS))
# dbRDA1 - 15%
# dbRDA2 - 0.4%

# RDA biplot 1 and 2
plot(dbMEM.rdaS, scaling = 2, main = "", type = "none", xlab = c("db-RDA-1"),
     ylab = c("db-RDA-2"), xlim = c(-40, 40), ylim = c(-10, 10))
# Add points
points(AA1[,1]*3, AA1[,2]*3, col = "black", bg = "red", pch = 21, cex = 1.2) 
points(AL1[,1]*3, AL1[,2]*3, col = "black", bg = "blue", pch = 21, cex = 1.2) 
points(AA2[,1]*3, AA2[,2]*3, col = "black", bg = "green", pch = 21, cex = 1.2) 
points(AL2[,1]*3, AL2[,2]*3, col = "black", bg = "purple", pch = 21, cex = 1.2) 
arrows(0,0, scores(dbMEM.rdaS, display = "bp", choices = 1, scaling = 1)*70,
       scores(dbMEM.rdaS, display = "bp", choices = 2, scaling = 1)*70,
       col = "black", length = 0.1)
text(scores(dbMEM.rdaS, display = "bp", choices = 1, scaling = 1)*75,
     scores(dbMEM.rdaS, display = "bp", choices = 2, scaling = 1)*75,
     labels = c("MEM12", "Overhang", "Outcrop position"), col = "black", cex = 0.8, pos = 1)


# Is there an effect xy #### 
xy.rda <- rda(gen.imp ~ x + y, data = metadata)
signif.xy.rda <- anova.cca(xy.rda)
signif.xy.rda  # not significant

# dbMEM with normal rda ####
MEM.rda <- rda(gen.imp ~ ., data = Y, scale = T)
vif.cca(MEM.rda) # When VIF is greater than around 10, this is problematic.

RsquareAdj(MEM.rda) # 10%
summary(eigenvals(MEM.rda, model = "constrained")) # RDA1 60%, RDA2 4%
summary(eigenvals(MEM.rda, model = "unconstrained")) # PC1 41%, PC2 4.7%
summary(MEM.rda)

#Partitioning of correlations:
#Inertia Proportion
#Total           21764     1.0000
#Constrained      8022     0.3686
#Unconstrained   13742     0.6314

# Significance
signif.MEM.rda <- anova.cca(MEM.rda)
signif.MEM.rda  # significant

#signif.term.MEM.rda <- anova.cca(MEM.rda, by = 'terms')
#signif.term.MEM.rda # too hard for computer

MEM.rda0 <- rda(gen.imp ~ 1, Y)
Sel <- ordiR2step(MEM.rda0, scope = formula(MEM.rda), direction = 'both')
Sel$anova #only overhang and prop
ade4::s.value(coords[,1:2], dbMEM[,12]) # d = 5?

# MEM5, MEM1, MEM6, MEM12
# forward.sel(gen.imp, X)
forward.sel.par(gen.imp, dbMEM)
