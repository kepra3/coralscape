setwd("~/git/coralscape/results/")

struc.complex1 <- read.csv("sample_metadata_cur_kal_05m_20200214_decvis_02.ply.csv")
struc.complex2 <- read.csv("sample_metadata_cur_kal_10m_20200214_decvis_02.ply.csv")
struc.complex3 <- read.csv("sample_metadata_cur_kal_20m_20200214_decvis_02.ply.csv")
struc.complex4 <- read.csv("sample_metadata_cur_sna_10m_20200303_decvis_02.ply.csv")
struc.complex5 <- read.csv("sample_metadata_cur_sna_20m_20200303_decvis_02.ply.csv")
struc.complex6 <- read.csv("sample_metadata_cur_sna_05m_20200303_decvis_02.ply.csv")
struc.complex <- rbind(struc.complex1,struc.complex2,struc.complex3,struc.complex4,struc.complex5,struc.complex6)
struc.complex$colony_rugosity <- gsub("[", "", struc.complex$colony_rugosity, fixed = TRUE)
struc.complex$colony_rugosity <- gsub("]", "", struc.complex$colony_rugosity, fixed = TRUE)
write.csv(struc.complex, file = "~/git/coralscape/results/sample_metadata.csv", quote = FALSE, row.names = FALSE)

annotations1 <- read.csv("scaled_annotations_cur_kal_05m_20200214_decvis_02.ply.csv")
annotations2 <- read.csv("scaled_annotations_cur_kal_10m_20200214_decvis_02.ply.csv")
annotations3 <- read.csv("scaled_annotations_cur_kal_20m_20200214_decvis_02.ply.csv")
annotations4 <- read.csv("scaled_annotations_cur_sna_10m_20200303_decvis_02.ply.csv")
annotations5 <- read.csv("scaled_annotations_cur_sna_20m_20200303_decvis_02.ply.csv")
annotations6 <- read.csv("scaled_annotations_cur_sna_05m_20200303_decvis_02.ply.csv")
annotations7 <- read.csv("scaled_annotations_cur_cas_05m_20201212_decvis_02.ply.csv")
annotations <- rbind(annotations1,annotations2, annotations3, annotations4, annotations5, annotations6, annotations7)
write.csv(annotations, file = "~/git/coralscape/results/annotations.csv", quote = FALSE, row.names = FALSE)

struc.complex4 <- read.csv("sample_metadata_cur_sna_10m_20200303_decvis_02.ply.csv")
struc.complex4_1 <- read.csv("sample_metadata_cur_sna_10m_20200303_decvis_02.ply_thetacompare.csv")
boxplot(struc.complex4$xz)
boxplot(struc.complex4_1$xz)
boxplot(struc.complex4$yz)
boxplot(struc.complex4_1$yz) # does normalise it
# might still be having the same problem...

boxplot(struc.complex3$overhang_prop)

boxplot(struc.complex5$overhang_prop)
boxplot(struc.complex5$colony_theta_xz)
boxplot(struc.complex5$colony_psi_yz)
struc.complex5$colony_rugosity <- gsub("[", "", struc.complex5$colony_rugosity, fixed = TRUE)
struc.complex5$colony_rugosity <- gsub("]", "", struc.complex5$colony_rugosity, fixed = TRUE)
struc.complex <- rbind(struc.complex3,struc.complex5)
write.csv(struc.complex, file = "~/git/coralscape/results/sample_metadata_X.csv", quote = FALSE, row.names = FALSE)

struc.complex_WP20 <- read.csv("struc_complex_results_WP20.txt", sep = "\t")
struc.complex_results <- read.csv("struc_complex_results.txt", sep = "\t")
struc.complex.X <- rbind(struc.complex_results, struc.complex_WP20)
write.csv(struc.complex.X, file = "struc_complex_results_X.csv", quote = FALSE, row.names = FALSE)

# Cluster file ####
all.aga <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/all-aga_1d_nc_20_6.csv")
hu <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/hu_1d_nc_20_4.csv")

order.all.aga <-  all.aga[order(all.aga[,1]),]
order.hu <- hu[order(hu[,1]),]
order.all.aga$hu.clust <- NA
order.all.aga$hu.clust[order.all.aga$Individual %in% order.hu$Individual] <- order.hu$Clusters

for (i in 1:length(order.all.aga[,1])) {
  if (is.na(order.all.aga$hu.clust[i])) {
  print('not hu')
  } else if (order.all.aga$hu.clust[i] == "Clust3") {
    order.all.aga$Clusters[i] = 7
    print('changing cluster to 7')
  } else {
    print('another cluster')
  }
}
write.csv(order.all.aga[,-20], file = "clusters_updated.csv", quote = FALSE, row.names = FALSE)
