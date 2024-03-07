###CHAPTER 2 GRAPHS AND FIGURES###
library(tidyverse)
library(tidylog)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(multcomp)
library(multcompView)

##FD metrics accross aridity and grazing####
#import FD results
FD_results <- read.csv("Functional trait data\\results\\FD_results_4Mar2024.csv", row.names = 1)
FD_results$ID <- as.factor(FD_results$ID)
#add grazing and aridity
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(c(ID, GRAZ, ARIDITY.v3))
siteinfo$ID <-as.factor(siteinfo$ID)
#do the join
long_FD_results <- FD_results |> 
  inner_join(siteinfo, by = "ID") |> 
  select(!c(qual.FRic, RaoQ, nsp)) |> 
  pivot_longer(cols = c(FRic, FEve, FDiv), names_to = "FD_metrics", values_to = "value" )
long_FD_results$GRAZ <- as.factor(long_FD_results$GRAZ)


FD_arid_grad <- ggplot(long_FD_results, aes(y = value, x = ARIDITY.v3)) +
  geom_point() +
  theme_classic() +
  xlab("Aridity") +
  facet_wrap(~ FD_metrics, scales = "free_y")

FD_graz_grad <- ggplot(long_FD_results, aes(y = value, x = GRAZ, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  theme_classic() +
  scale_fill_manual(values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  xlab("") +
  scale_x_discrete(labels = c("Ungrazed", "Low grazing", "Medium grazing", "High grazing")) +
  facet_wrap(~ FD_metrics, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank(), legend.position = "none")

FD_grad_plots <- ggarrange(FD_arid_grad, FD_graz_grad, ncol = 1, nrow= 2, labels = c("a", "b"))
ggsave("FD_grad_plots.png", FD_grad_plots, height = 1600, width = 1300, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")



##Nint ~ FD metrics####
#import nint results
#Add the NIntc
nint_result <- 
  read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1) |> 
  filter(ID %in% c(FD_results$ID))
nint_result$ID <- as.factor(nint_result$ID)


##summarise the nintc by plot
nint_sum <- nint_result |> 
  select(country, site_ID, ID, graz, aridity, NIntc_richness, NIntc_cover, NInta_richness, NInta_cover) |> 
  group_by(ID) |> 
  #calculate mean and sd of NIntc richness
  mutate(mean_NIntc_rich = mean(NIntc_richness, na.rm = TRUE), 
         se_NIntc_rich = sd(NIntc_richness, na.rm = TRUE)/sqrt(sum(!is.na(NIntc_richness)))) |> 
  #calculate the mean and sd of NIntc cover
  mutate(mean_NIntc_cov = mean(NIntc_cover, na.rm = TRUE), 
         se_NIntc_cov = sd(NIntc_cover, na.rm = T)/sqrt(sum(!is.na(NIntc_cover)))) |> 
  #calculate mean and sd of NInta richness
  mutate(mean_NInta_rich = mean(NInta_richness, na.rm = T), 
         se_NInta_rich = sd(NInta_richness, na.rm = T)/sqrt(sum(!is.na(NInta_richness)))) |> 
  #calculate the mean and sd of NInta cover
  mutate(mean_NInta_cov = mean(NInta_cover, na.rm = T), 
         se_NInta_cov = sd(NInta_cover, na.rm = T)/sqrt(sum(!is.na(NInta_cover)))) |>
  ungroup() |> 
  select(!c(NIntc_richness, NIntc_cover, NInta_richness, NInta_cover)) |> 
  distinct() |> #remove duplicate rows, only need one eman per plot
  left_join(FD_results, by = "ID") |>  #join to the FD_results
  filter(!is.na(FEve))


##NIntc richness against the FD metrics##
rich_fdiv <- ggplot(nint_sum, aes(x = FDiv, y = mean_NIntc_rich)) +
  geom_point() +
  ylab("") +
  geom_errorbar(aes(ymin = mean_NIntc_rich - se_NIntc_rich, ymax =  mean_NIntc_rich + se_NIntc_rich), width = 0.01) +
  theme_classic()

rich_frich <- ggplot(nint_sum, aes(x = FRic, y = mean_NIntc_rich)) +
  geom_point() +
  ylab("mean NIntc richness") +
  geom_errorbar(aes(ymin = mean_NIntc_rich - se_NIntc_rich, ymax =  mean_NIntc_rich + se_NIntc_rich), 
  width = 0.01) +
  theme_classic()

rich_feve <- ggplot(nint_sum, aes(x = FEve, y = mean_NIntc_rich)) +
  geom_point() +
  ylab("") +
  geom_errorbar(aes(ymin = mean_NIntc_rich - se_NIntc_rich, ymax =  mean_NIntc_rich + se_NIntc_rich), 
                width = 0.01) +
  theme_classic()

nintc_rich_FD <- ggarrange(rich_frich, rich_feve, rich_fdiv, ncol = 3, nrow = 1, labels = c("a", "b", "c"))
ggsave("NIntc_FD_scatterplots.png", nintc_rich_FD, height = 800, width = 1700, units = "px", 
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


###Functional distance ~ GRAZ, ARIDITY, microsite affinty####
#import data
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp.csv", row.names = 1)
twosp_dist$GRAZ <- as.factor(twosp_dist$GRAZ)
twosp_dist$SITE_ID <- as.factor(twosp_dist$SITE_ID)
twosp_dist$grouping <- as.factor(twosp_dist$grouping)

#the model
twosp_dist$arid_sq <- (twosp_dist$ARIDITY.v3)^2
twosp2 <- glmmTMB(euclidean_dist ~ grouping*GRAZ +grouping*ARIDITY.v3 + 
                    ARIDITY.v3*GRAZ + grouping*arid_sq + arid_sq*GRAZ + (1|SITE_ID), data = twosp_dist)
Anova(twosp2, type = 2) #type 2 is seq SS
summary(twosp2)

cld(glht(model = twosp2, mcp(GRAZ = "Tukey")))
emmeans(twosp2, specs = "GRAZ")

##distance~grouping, colour by graz
#we want to know which levels of graz are significantly different from each other, within gspecific groupings
#concateante graz and grouping to create treatemnt ID
twosp_dist <- twosp_dist |> 
  mutate(treatmentID = str_c(grouping, GRAZ, sep = "_"))
twosp_dist$treatmentID <- as.factor(twosp_dist$treatmentID)


treatment_mod <- lm(euclidean_dist ~ treatmentID, data = twosp_dist)
Anova(treatment_mod)
letters <- cld(glht(model = treatment_mod, mcp(treatmentID = "Tukey")))
#look at letters and make a vector to of them to add to the graph
#put them in the same order as the boxplots on the graph
letters_forplot <- data.frame(letters = c("a", "bc", "cde", "ef", "ag", "bc", "ef", "ce", "af", "b", "bd", "efg"), 
                              xcoord = c(rep("bare", 4), rep("dominant", 4), rep("both", 4)),
                              ycoord = c(rep(9, 12)))

dist_bygrouping_graz <- ggplot(twosp_dist, aes(x = grouping, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(labels = c("ungrazed", "low", "medium", "high"), 
                    values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  labs(y = "Euclidean distance between dominant and target species", 
       x = "Microsite affinity of target species", 
       fill = "") +
  scale_x_discrete(labels = c("bare", "dominant", "both")) +
  geom_text(aes(x = GRAZ, y = 9, ),label = c(letters_forplot$letters)) + ##fix this
  theme_classic() +
  theme(legend.position = "bottom") 

ggsave("dist_grouping_graz.png", dist_bygrouping_graz, height = 1500, width = 1200, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


##distance~aridity
loess_color <- brewer.pal(8, "Dark2")[7]
dist_aridity <- ggplot(twosp_dist, aes(x = ARIDITY.v3, y = euclidean_dist)) +
  geom_jitter(width = 0.08, height = 0.2, alpha = 0.6) +
  geom_smooth(method = "lm", color = loess_color, fill = loess_color) +
  facet_wrap(~grouping, ncol = 1,
             labeller = as_labeller(c("nurse_bare_only" = "bare", "nurse_both" = "both", "nurse_nurse_only" = "nurse")))+
  labs(y = "Euclidean distance between dominant and target species", x = "Aridity") +
  theme_classic()

ggsave("dist_aridity.png", dist_aridity, height = 1400, width = 900, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")
