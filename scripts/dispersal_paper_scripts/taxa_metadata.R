# Organising taxa metadata
# Created 25/5/22

# Functions
find_taxa_4clones <- function(pop, clonegroups) {
  clonegroups$Taxa <- pop$Taxa[match(clonegroups$Sample, pop$Individual)]
  clonegroup.list <- clonegroups[is.na(clonegroups$Taxa),]$Groups
  for (i in clonegroup.list) {
    mlg.taxa <- clonegroups[clonegroups$Groups == i & !is.na(clonegroups$Taxa),4]
    clonegroups$Taxa[clonegroups$Groups == i] = mlg.taxa
  }
  return(clonegroups)
}

# Import datasets
setwd("~/git/coralscape/results/")
struc.complex <- read.csv("struc_complex_results_envonly.txt", sep = "\t")
taxa <- read.csv("all-aga_1d_wc_20_4_assigned.csv")
coordinates <- read.csv("~/git/coralscape/results/annotations.csv")

pop_ac_1div_nc <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/ac_1div_nc_20_4.csv")
pop_hu_1div_nc <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/hu_1div_nc_20_4.csv")
pop_lm_1div_nc <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/lm_1div_nc_20_2.csv")
ac.clonegroups <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/Pyscripts/results/clones/clone_groups_ac.csv")
hu.clonegroups <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/Pyscripts/results/clones/clone_groups_hu.csv")
lm.clonegroups <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/Pyscripts/results/clones/clone_groups_lm.csv")

pop_ac_1div_nc$Taxa[is.na(pop_ac_1div_nc$Taxa)]
pop_hu_1div_nc$Taxa[is.na(pop_hu_1div_nc$Taxa)] <- "admixed"
pop_lm_1div_nc$Taxa[is.na(pop_lm_1div_nc$Taxa)]

# Remove 529 as a hybrid between lm and ac
# outlier individuals already removed
lm.clonegroups[lm.clonegroups$Sample=="KP0529_LM_WP20",] # group 93
lm.clonegroups[lm.clonegroups$Sample=="KP0631_LM_SB20",] # group 91
lm.clonegroups[lm.clonegroups$Sample=="KP1093_LM_CA20",] # group 92
lm.clonegroups[lm.clonegroups$Sample=="KP0670_LM_SB20",] # group 90
lm.clonegroups <- lm.clonegroups[-c(35, 54, 87, 62),]

ac.clonegroups <- find_taxa_4clones(pop_ac_1div_nc, ac.clonegroups)
hu.clonegroups <- find_taxa_4clones(pop_hu_1div_nc, hu.clonegroups)
lm.clonegroups <- find_taxa_4clones(pop_lm_1div_nc, lm.clonegroups)

ac.clonegroups$Groups <- paste0("aa.", ac.clonegroups$Groups)
ac.clonegroups$Taxa[ac.clonegroups$Taxa == "Clust1"] <- "AA1"
ac.clonegroups$Taxa[ac.clonegroups$Taxa == "Clust2"] <- "AA2"
hu.clonegroups$Groups <- paste0("ah.", hu.clonegroups$Groups)
hu.clonegroups$Taxa[hu.clonegroups$Taxa == "Clust1"] <- "AH1"
hu.clonegroups$Taxa[hu.clonegroups$Taxa == "Clust2"] <- "AH2"
hu.clonegroups$Taxa[hu.clonegroups$Taxa == "Clust3"] <- "AH3"
lm.clonegroups$Groups <- paste0("al.", lm.clonegroups$Groups)
lm.clonegroups$Taxa[lm.clonegroups$Taxa == "Clust1"] <- "AL1"
lm.clonegroups$Taxa[lm.clonegroups$Taxa == "Clust2"] <- "AL2"

all.clonegroups <- rbind(ac.clonegroups, hu.clonegroups, lm.clonegroups)

taxa$Taxa <- all.clonegroups$Taxa[match(taxa$Individual, all.clonegroups$Sample)]
taxa$Taxa[is.na(taxa$Taxa)] <- "admixed"
taxa$Groups <- all.clonegroups$Groups[match(taxa$Individual, all.clonegroups$Sample)]

taxa$Groups[taxa$Individual=="KP0529_LM_WP20"] <- 93 # group 93
taxa$Groups[taxa$Individual=="KP0631_LM_SB20"] <- 91 # group 91
taxa$Groups[taxa$Individual=="KP1093_LM_CA20"] <- 92 # group 92
taxa$Groups[taxa$Individual=="KP0670_LM_SB20"] <- 90# group 90

write.csv(taxa, "taxa_metadata.csv", quote = FALSE, row.names = FALSE)
