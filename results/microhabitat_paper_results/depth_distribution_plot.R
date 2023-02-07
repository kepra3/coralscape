# Abundance and distribution of samples

ac_3b_nc_20_AA1_all_ <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/3b - Within & between plots/spatial-agaricia/data/ac_3b_nc_20_AA1_all_.csv", row.names = 1)


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
ggsave("depth_distribution.pdf", height = 4, width = 17, units = "cm")