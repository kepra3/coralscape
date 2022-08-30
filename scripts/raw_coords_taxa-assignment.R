clean_df <- function(data) {
  data <- data %>% separate(X, into = c("Sample", "Sp", "Site", "extra"), sep = "_", remove = TRUE) %>% 
    unite(Individual, Sample, Sp, Site)
  return(data)
}

# Improt data
df <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/Photogrammetry/CloudCompare/WP20/cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt", header=FALSE)
taxa.ac <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/ac_1div_nc_20_4.csv")
taxa.hu <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/hu_1div_nc_20_4.csv")
taxa.lm <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/lm_1div_nc-wnr_20_2.csv")

# Clean df
colnames(df) <- c("X", "x", "y", "z")
df <- clean_df(df)
df <- df[is.na(df$extra),]
df <- df %>% dplyr::select(Individual, x, y, z)

# Set taxa names
taxa.ac$Taxa[taxa.ac$Taxa == "Clust2"] = "AA2"
taxa.ac$Taxa[taxa.ac$Taxa == "Clust1"]  = "AA1"
taxa.hu$Taxa[taxa.hu$Taxa == "Clust1"] = "AH1"
taxa.hu$Taxa[taxa.hu$Taxa == "Clust2"]  = "AH2"
taxa.hu$Taxa[taxa.hu$Taxa == "Clust3"] = "AH3"
taxa.lm$Taxa[taxa.lm$Taxa == "Clust1"] = "AL1"
taxa.lm$Taxa[taxa.lm$Taxa == "Clust2"]  = "AL2"

ac.df <- df
hu.df <- df
lm.df <- df
ac.df$Taxa <- taxa.ac$Taxa[match(df$Individual, taxa.ac$Individual)]
hu.df$Taxa <- taxa.hu$Taxa[match(df$Individual, taxa.hu$Individual)]
lm.df$Taxa <- taxa.lm$Taxa[match(df$Individual, taxa.lm$Individual)]

ac.df <- na.omit(ac.df)
hu.df <- na.omit(hu.df)
lm.df <- na.omit(lm.df)

colours <- c("#274e13","#8fce00", "#b06100", "#ffa500", "#ff4e00", "#6f0641", "#7223f0")

ac.df$colours[ac.df$Taxa == "AA1"] <- colours[1]
ac.df$colours[ac.df$Taxa == "AA2"] <- colours[2]
hu.df$colours[hu.df$Taxa == "AH1"] <- colours[3]
hu.df$colours[hu.df$Taxa == "AH2"] <- colours[4]
hu.df$colours[hu.df$Taxa == "AH3"] <- colours[5]
lm.df$colours[lm.df$Taxa == "AL1"] <- colours[6]
lm.df$colours[lm.df$Taxa == "AL2"] <- colours[7]

all.df <- rbind(ac.df, hu.df, lm.df)
colnames(all.df) <- c("label", "x", "y", "z", "taxa", "color")
all.df <- all.df[, c(2:4, 1, 6)]
write.csv(all.df, file = "~/Dropbox/agaricia_project_2019/shalo_ag/cur_kal_20m_20201214_agaricia-annotations-KP.csv", quote = FALSE, row.names = FALSE)
