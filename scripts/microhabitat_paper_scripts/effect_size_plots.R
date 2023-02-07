# Making effect size plots

# Packages
library(ggplot2)

# Working directory
setwd("~/git/coralscape/results")

effect_size_plot <- function(microhabitat_variable) {
  five <- read.csv(paste0("univariate_emmeans/", microhabitat_variable, "5.emmeans.txt"), row.names = 1)
  ten <- read.csv(paste0("univariate_emmeans/", microhabitat_variable, "10.emmeans.txt"), row.names = 1)
  twenty <- read.csv(paste0("univariate_emmeans/", microhabitat_variable, "20.emmeans.txt"), row.names = 1)
  five$Depth <- "5"
  ten$Depth <- "10"
  twenty$Depth <- "20"
  all <- rbind(five, ten, twenty)
  all$Depth <- factor(all$Depth, levels = c("5", "10", "20"))
  # Swapping estimate diff so that deeper colony is always first
  all$estimate[all$X1 == "AA1 - AL1"] = -(all$estimate[all$X1 == "AA1 - AL1"])
  all$X1[all$X1 == "AA1 - AL1"] = "AL1 - AA1"
  
  all$estimate[all$X1 == "AA1 - AL2"] = -(all$estimate[all$X1 == "AA1 - AL2"])
  all$X1[all$X1 == "AA1 - AL2"] = "AL2 - AA1"
  
  all$estimate[all$X1 == "AA2 - AL1"] = -(all$estimate[all$X1 == "AA2 - AL1"])
  all$X1[all$X1 == "AA2 - AL1"] = "AL1 - AA2"
  
  all$estimate[all$X1 == "AA2 - AL2"] = -(all$estimate[all$X1 == "AA2 - AL2"])
  all$X1[all$X1 == "AA2 - AL2"] = "AL2 - AA2"
  
  all$estimate[all$X1 == "AH1 - AL1"] = -(all$estimate[all$X1 == "AH1 - AL1"])
  all$X1[all$X1 == "AH1 - AL1"] = "AL1 - AH1"
  
  all$estimate[all$X1 == "AH1 - AL2"] = -(all$estimate[all$X1 == "AH1 - AL2"])
  all$X1[all$X1 == "AH1 - AL2"] = "AL2 - AH1"
  
  all$estimate[all$X1 == "AH3 - AL2"] = -(all$estimate[all$X1 == "AH3 - AL2"])
  all$X1[all$X1 == "AH3 - AL2"] = "AL2 - AH3"
  # Order of comparisons
  all$X1 <- factor(all$X1, levels = c("AA1 - AA2", # deep to shallow
                                      "AH1 - AH2", "AH1 - AH3", "AH2 - AH3",# deep to shallow
                                      "AL1 - AL2", # deep to shallow
                                      "AA1 - AH1","AA2 - AH1", # deep to shallow
                                      "AA2 - AH2", "AA2 - AH3", # deep to shallow
                                      # swapped all these, now deep to shallow
                                      "AL1 - AA1", "AL2 - AA1", 
                                      "AL1 - AA2", "AL2 - AA2",
                                      "AL1 - AH1", "AL2 - AH1", "AL2 - AH3"))
  if (microhabitat_variable == "overhang_prop") {
    for (se_est in 1:length(all$SE)) {
      if (all$SE[se_est] > 10) {
        all$SE[se_est] = NA
        all$estimate[se_est] = NA
      }
    }
  }
 
  # Plot results
  plot <- ggplot(data = all, aes(x = X1, y = estimate)) + geom_point() + 
    geom_errorbar(aes(y = estimate, ymin=estimate-SE, ymax=estimate+SE)) + 
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    scale_x_discrete(position = "top") +
    facet_wrap(~Depth, nrow = 3, strip.position = "left") +
    theme_bw() +
    xlab(paste0("Pairwise comparisons - ", microhabitat_variable)) +
    ylab("Effect size + SE") +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.text.x.top = element_text(angle = 90))
  return(plot)
}

ground_elevation_plot <- effect_size_plot("ground_elevation")
ground_elevation_plot
ggsave("ground_elevation_effect_size.pdf", width = 10, height = 12, units = "cm")
env_rugosity_plot <- effect_size_plot("environment_rugosity")
env_rugosity_plot
ggsave("environment_rugosity_effect_size.pdf", width = 10, height = 12, units = "cm")
outcrop_proportion_plot <- effect_size_plot("outcrop_prop")
outcrop_proportion_plot
ggsave("outcrop_prop_effect_size.pdf", width = 10, height = 12, units = "cm")
overhang_proportion_plot <- effect_size_plot("overhang_prop")
overhang_proportion_plot
ggsave("overhang_prop_effect_size.pdf", width = 10, height = 12, units = "cm")
local_height_plot <- effect_size_plot("z_scale")
local_height_plot
ggsave("z_scale_effect_size.pdf", width = 10, height = 12, units = "cm")
