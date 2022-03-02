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
annotations <- rbind(annotations1,annotations2, annotations3, annotations4, annotations5, annotations6)
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

