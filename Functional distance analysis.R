##Functional distance analysis
library(tidyverse)
library(tidylog)
library(ggplot2)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(MuMIn)

#From the filled trait data for plotspecific species
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv", row.names = 1)

##We need to standardise each trait value to mean 0 and variance of 1 with Z = (x-mean)/sd
##We need only one trait value per species

#Get the mean of each trait for each sp
FT_mean <- FT |> 
  group_by(taxon,trait) |> 
  summarise(mean_value = mean(value))

#Get it into wide format
FT_wide <- FT_mean |>
  pivot_wider(names_from = trait, values_from = mean_value) |> 
  column_to_rownames(var = "taxon") |> 
  filter(!is.na(MaxH),
         !is.na(MaxLS), 
         !is.na(MeanLA),
         !is.na(MeanLDMC),
         !is.na(MeanLL),
         !is.na(MeanSLA), 
         !is.na(percentC),
         !is.na(percentN)) |> 
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN))



#standardise trait values
std_FT_wide <- FT_wide

traitlist <- c("MaxH","MaxLS","MeanLA","MeanLDMC","MeanLL","MeanSLA","C_N_ratio")

for(t in 1:length(traitlist)) {
  
  #get the grand mean and sd for each trait
  m <- FT_wide |> 
    summarise(m = mean(FT_wide[ , which(colnames(FT_wide) == traitlist[t])]))
  
  sd <- FT_wide |> 
    summarise(m = sd(FT_wide[ , which(colnames(FT_wide) == traitlist[t])]))
  
  for (i in 1:nrow(FT_wide)) {
    
    raw_value <- FT_wide[i , which(colnames(FT_wide) == traitlist[t])]
    #replace the raw value with the standardised value using m and sd
    std_FT_wide[i , which(colnames(std_FT_wide) == traitlist[t])] <- (raw_value - m)/sd
    
  }
}


#Get the euclidean distance between each pair of species
distmat <- as.matrix(dist(std_FT_wide, method = "euclidean"))

##Get the species that are nurses, growing with nurses, or in the bare microsite in each country####
#We require raw country data
data_files <- list.files("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3\\",
                          data_files[i])))
}

countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")

##LOOP starts here
l = 1
for (t in 1:length(countrynames)) {
  cou <- get(countrynames[t])
  
  IDlist <- unique(cou$ID)
  
    for(i in 1:length(IDlist)) {
    plot <- cou |> 
      filter(ID == IDlist[i])
    
    replist <- unique(plot$Number.of.replicate)
  
      for(r in 1:length(replist)) {
        one_rep <- plot |> 
          filter(Number.of.replicate == replist[r])
        
        NURSE <- one_rep |> 
          filter(Microsite == 2) |> 
          distinct(ID_Microsite)
        
        bare_sp <- one_rep |> #sp in bare microsite
          filter(Microsite == 1) |> 
          select(Species.within.quadrat)
        
        nurse_sp <- one_rep |> #sp in nurse microsite
          filter(Microsite == 2) |> 
          select(Species.within.quadrat)
        
        #get the species that are in both nurse and bare microsites
        both_sp <- c(bare_sp[match(bare_sp$Species.within.quadrat, nurse_sp$Species.within.quadrat) , ])
        both_sp <- both_sp[!is.na(both_sp)]
        #if there are no species in both microsites, assign NA
        if(length(both_sp) == 0) {
          both_sp <- NA
        }
        
        #only if both_sp is not NA:
        if(FALSE %in% is.na(both_sp)) {
        
        #species ONLY in nurse microsite:
        nurse_only <- nurse_sp[-which(nurse_sp$Species.within.quadrat %in% both_sp), ]
        #if there are no species, assign NA
        if(length(nurse_only) == 0) {
          nurse_only <- NA
        }
        
        #species ONLY in bare microsites
        bare_only <- bare_sp[-which(bare_sp$Species.within.quadrat %in% both_sp) , ]
        #if there are no species, assign NA
        if(length(bare_only) == 0) {
          bare_only <- NA
        }
        
        #if both sp is NA:
        }else {
          nurse_only <- nurse_sp$Species.within.quadrat
          bare_only <- bare_sp$Species.within.quadrat
        }
        
        #only for the first run of the loop:
        if(l == 1){
        #put species and their classifications in the dataframe
        nurse_identity_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = NURSE$ID_Microsite, position = "nurse_species")
        nurse_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = nurse_only, position = "nurse_only")
        bare_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = bare_only, position = "bare_only")
        both_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = both_sp, position = "both")
        
        sp_positions <- rbind(nurse_identity_df, nurse_df, bare_df, both_df)
        
        } else { #for subsequent runs we rbind to the df we made whn l was 1
          nurse_identity_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = NURSE$ID_Microsite, position = "nurse_species")
          nurse_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = nurse_only, position = "nurse_only")
          bare_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = bare_only, position = "bare_only")
          both_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = both_sp, position = "both")
          
          sp_positions <- rbind(sp_positions, nurse_identity_df, nurse_df, bare_df, both_df)
        }
        
        l = l+1
      }#close loop through replicates
  }#close loop through ID's
}#close loop through countries

sp_positions <- as.data.frame(sp_positions)
#save the table 
write.csv(sp_positions, "Functional trait data\\results\\sp_positions.csv")


###Retreive distance between pairs from distance matrix####
#import sp_positions
sp_positions <- read.csv("Functional trait data\\results\\sp_positions.csv", row.names = 1) 

IDlist <- unique(sp_positions$ID)

pair_dist <- data.frame(ID = character(), replicate = character(), 
                        D_bare_only = numeric(), D_nurse_only = numeric(), D_both = numeric())

l = 1
for(i in 1:length(IDlist)) {
  plot <- sp_positions |> 
    filter(ID == IDlist[i])
  
  replist <- unique(plot$replicate) 
  
  for(r in 1:length(replist)) {
    #Get the names of nurse, bare only, nurse only and both species
    NURSE <- sp_positions |> 
      filter(ID == IDlist[i], replicate == replist[r], position == "nurse_species") |> 
      select(taxon)
    
    bare_only <- sp_positions |> 
      filter(ID == IDlist[i], replicate == replist[r], position == "bare_only") |> 
      select(taxon)
    
    nurse_only <- sp_positions |> 
      filter(ID == IDlist[i], replicate == replist[r], position == "nurse_only") |> 
      select(taxon)
    
    both <- sp_positions |> 
      filter(ID == IDlist[i], replicate == replist[r], position == "both") |> 
      select(taxon)
    
    #Get mean distance between the NURSE and bare_only, nurse_only and both species
    #mean dist between nurse and sp growing only in bare microsite
    D_bare_only <- 
      mean(distmat[which(rownames(distmat) == NURSE$taxon), which(colnames(distmat) %in% c(bare_only$taxon))])
    #mean dist between nurse and sp growing only in nurse microsite
    D_nurse_only <- 
      mean(distmat[which(rownames(distmat) == NURSE$taxon), which(colnames(distmat) %in% c(nurse_only$taxon))])
    #mean distance between nurse and sp growing in both microsites
    D_both <- 
      mean(distmat[which(rownames(distmat) == NURSE$taxon), which(colnames(distmat) %in% c(both$taxon))])
    
    pair_dist[l, 1] <- IDlist[i]
    pair_dist[l, 2] <- replist[r]
    pair_dist[l, 3] <- D_bare_only
    pair_dist[l, 4] <- D_nurse_only
    pair_dist[l, 5] <- D_both
    
    l = l+1
  }
}
#The distance is NaN if nurse_only, bare_only or both is NA
write.csv(pair_dist, "Functional trait data\\results\\Functional_distances.csv")


##Now we need to add the aridity, graz of each plot, and the nintc of each rep
pair_dist <- read.csv("Functional trait data\\Results\\Functional_distances.csv", row.names = 1)
pair_dist$ID <- as.numeric(pair_dist$ID)
pair_dist$replicate <- as.numeric(pair_dist$replicate)

#import siteinfo
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(ID,COU,SITE,SITE_ID,PLOT,GRAZ, ARIDITY.v3) |> 
  distinct() |> 
  filter(!is.na(ID))
#do the join
pair_dist <- pair_dist |> 
  inner_join(siteinfo, by = "ID") |> 
  rename(replicate_no = replicate)
pair_dist$GRAZ <- as.factor(pair_dist$GRAZ)
pair_dist$SITE_ID <- as.factor(pair_dist$SITE_ID)


#Add the NIntc
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1) |> 
  select(ID, replicate_no, NIntc_richness, NIntc_cover, NInta_richness, NInta_cover)

pair_dist <- pair_dist |> 
  left_join(nint_result, by = c("ID", "replicate_no")) |> 
  mutate(Delta_nurse_bare = D_nurse_only - D_bare_only)

ggplot(pair_dist, aes(x = Delta_nurse_bare, y = NIntc_richness)) +
  geom_jitter(width = 1, height = 1)

ggplot(pair_dist, aes(x = D_nurse_only, y = NIntc_richness)) +
  geom_jitter(width = 1, height = 1)

ggplot(pair_dist, aes(x = D_bare_only, y = NIntc_richness)) +
  geom_jitter(width = 1, height = 1)

ggplot(pair_dist, aes(y = Delta_nurse_bare, x = ARIDITY.v3)) +
  geom_jitter(width = 1, height = 1)

ggplot(pair_dist, aes(y = Delta_nurse_bare, x = GRAZ)) +
 geom_boxplot()

##Let's make some models
#NIntc is bounded beween -1 and 1, so binomial family is appropriate
#However the function requires that the response be bounded between 0 and 1, so rescale NIntc
#x-min/max- min (here the formula is just already simplified)
pair_dist$NIntc_richness_binom <- (pair_dist$NIntc_richness + 1)/2
pair_dist$NIntc_cover_binom <- (pair_dist$NIntc_cover + 1)/2

#x-min/max- min
pair_dist$NInta_richness_binom <- (pair_dist$NInta_richness - (-1)) / (2 - (-1))
pair_dist$NInta_cover_binom <- (pair_dist$NInta_cover - (-1)) / (2 - (-1))

mod <- glmmTMB(NIntc_richness_binom ~ D_nurse_only + (1|SITE_ID), family = "binomial", data = pair_dist)
Anova(mod)
summary(mod)

mod2 <- glmmTMB(NIntc_richness_binom ~ Delta_nurse_bare + (1|SITE_ID), family = "binomial", data = pair_dist)
Anova(mod2)
summary(mod2)

mod3 <- glmmTMB(NIntc_richness_binom ~ D_bare_only + (1|SITE_ID), family = "binomial", data = pair_dist)
Anova(mod3)
summary(mod2)



#Does the functional distance depend on whether the distance is between nurse and nurse_only species, or nurse and bare_only species, or nurse and both species
#We need to change the pair_dist data to long format

long_pair_dist <- pair_dist |> 
  select(ID, SITE_ID, replicate_no, D_bare_only, D_nurse_only, D_both, GRAZ, 
         ARIDITY.v3, NIntc_richness_binom, NIntc_cover_binom, NInta_richness_binom, NInta_cover_binom) |> 
  pivot_longer(cols = c(D_bare_only, D_nurse_only, D_both), names_to = "grouping", values_to = "mean_Fdist")
long_pair_dist$grouping <- as.factor(long_pair_dist$grouping)

mod4 <- glmmTMB(mean_Fdist ~ grouping + GRAZ + ARIDITY.v3 + (1|SITE_ID), data = long_pair_dist)
Anova(mod4)
#Graz significant
summary(mod4)
emmeans(mod4, specs = c("GRAZ")) #?? does not agree with summary table

#without random effect
mod5 <- glmmTMB(mean_Fdist ~ grouping + GRAZ + ARIDITY.v3 , data = long_pair_dist)
Anova(mod5)
#now graz and aridity is significant
summary(mod5)
emmeans(mod5, specs = "GRAZ")


mod6 <- glmmTMB(mean_Fdist ~ grouping*GRAZ*ARIDITY.v3 + (1|SITE_ID), data = long_pair_dist)
Anova(mod6)
#Graz significant
summary(mod4)
emmeans(mod4, specs = c("GRAZ")) #?? does not agree with summary table

ggplot(long_pair_dist, aes(x = mean_Fdist)) +
  geom_histogram()

ggplot(long_pair_dist, aes(x = GRAZ, y = mean_Fdist)) +
  geom_boxplot() 

ggplot(long_pair_dist, aes(x = grouping, y = mean_Fdist, fill = GRAZ)) +
  geom_boxplot()

ggplot(long_pair_dist, aes(x = grouping, y = mean_Fdist)) +
  geom_boxplot()

ggplot(long_pair_dist, aes(x = ARIDITY.v3, y = mean_Fdist)) +
  geom_jitter(width = 0.1, height = 0.1)

#let's bin aridity
long_pair_dist <- long_pair_dist |> 
       mutate(arid_bin = cut(ARIDITY.v3, breaks=4))

ggplot(long_pair_dist, aes(x = arid_bin, y = mean_Fdist)) +
  geom_boxplot()


###Another way: retreive distances between the nurse and each species also growing in the rep####
#then classify that distance as D_bare_only, D_nurse_only or D_both
#This distance is only between 2 sp, not between the nurse and many sp as in pair_dist

#import sp_positions
sp_positions <- read.csv("Functional trait data\\results\\sp_positions.csv", row.names = 1) 

IDlist <- unique(sp_positions$ID)

twosp_dist <- cbind("ID" = NA, "replicate" = NA, 
                    "euclidean_dist" = NA, "grouping" = NA)

for(i in 1:length(IDlist)) {
  plot <- sp_positions |> 
    filter(ID == IDlist[i])
  
  replist <- unique(plot$replicate) 
  
  for(r in 1:length(replist)) {
    #Get the names of nurse, bare only, nurse only and both species
    NURSE <- sp_positions |> 
      filter(ID == IDlist[i], replicate == replist[r], position == "nurse_species") |> 
      select(taxon)
    NURSE <- NURSE$taxon
    
    bare_only <- sp_positions |> 
      filter(ID == IDlist[i], replicate == replist[r], position == "bare_only") |> 
      select(taxon)
    bare_only <- c(bare_only$taxon)
    
    nurse_only <- sp_positions |> 
      filter(ID == IDlist[i], replicate == replist[r], position == "nurse_only") |> 
      select(taxon)
    nurse_only <- c(nurse_only$taxon)
    
    both <- sp_positions |> 
      filter(ID == IDlist[i], replicate == replist[r], position == "both") |> 
      select(taxon)
    both <- c(both$taxon)
    
    ###
    for(b in 1:length(bare_only)) {
      D_bare_only <- 
      mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == bare_only[b])])
      
      if(b == 1){
        twosp_bare_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                 "euclidean_dist" = D_bare_only, "grouping" = "nurse_bare_only")
      } else {
        temp_twosp_bare_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                      "euclidean_dist" = D_bare_only, "grouping" = "nurse_bare_only")
        twosp_bare_only <- rbind(twosp_bare_only, temp_twosp_bare_only)
      }}
      
    ###
    for(n in 1:length(nurse_only)) {
      D_nurse_only <- 
        mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == nurse_only[n])])
      
      if(n == 1){
        twosp_nurse_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                 "euclidean_dist" = D_nurse_only, "grouping" = "nurse_nurse_only")
      } else {
        temp_twosp_nurse_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                      "euclidean_dist" = D_nurse_only, "grouping" = "nurse_nurse_only")
        twosp_nurse_only <- rbind(twosp_nurse_only, temp_twosp_nurse_only)
      }}
      
    ###
    for(z in 1:length(both)) {
      D_both <- 
        mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == both[z])])
      
      if(z == 1){
        twosp_both <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                  "euclidean_dist" = D_both, "grouping" = "nurse_both")
      } else {
        temp_twosp_both <- twosp_both <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                                     "euclidean_dist" = D_both, "grouping" = "nurse_both")
        twosp_both <- rbind(twosp_both, temp_twosp_both)
      }}
      
      twosp_dist <- rbind(twosp_dist, twosp_bare_only,twosp_nurse_only, twosp_both)
  }
}
#The distance is NaN if nurse_only, bare_only or both is NA

#Add ardidty and graz
twosp_dist <- as.data.frame(twosp_dist)
twosp_dist$ID <- as.numeric(twosp_dist$ID)
#import siteinfo
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(ID,COU,SITE,SITE_ID,PLOT,GRAZ, ARIDITY.v3) |> 
  distinct() |> 
  filter(!is.na(ID))

#do the join
twosp_dist <- twosp_dist |> 
  filter(!is.na(euclidean_dist), !euclidean_dist == "NaN") |> 
  inner_join(siteinfo, by = "ID") 

write.csv(twosp_dist, "Functional trait data\\results\\Functional_distances_between_2sp.csv")


###Make some models
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp.csv", row.names = 1)
twosp_dist$GRAZ <- as.factor(twosp_dist$GRAZ)
twosp_dist$SITE_ID <- as.factor(twosp_dist$SITE_ID)
twosp_dist$grouping <- as.factor(twosp_dist$grouping)

ggplot(twosp_dist, aes(x = euclidean_dist)) +
  geom_histogram()

twosp1 <- glmmTMB(euclidean_dist ~ grouping*GRAZ*ARIDITY.v3 + (1|SITE_ID), data = twosp_dist)
Anova(twosp1)
summary(twosp1)

#get whole model p
#null model:
twosp_null <- glmmTMB(euclidean_dist ~ 1 + (1|SITE_ID), data = twosp_dist)
anova(twosp_null, twosp1) #whole model p = < 2.2e-16


#get R squared
r.squaredGLMM(twosp1) #take the theoretical

#check residuals
twosp1_simres <- simulateResiduals(twosp1)
plot(twosp1_simres)

emmeans(twosp1, specs = "grouping")
emmeans(twosp1, specs = "GRAZ")


#model including the aridity squared 
twosp_dist$arid_sq <- (twosp_dist$ARIDITY.v3)^2
twosp2 <- glmmTMB(euclidean_dist ~ grouping*GRAZ +grouping*ARIDITY.v3 + 
                    ARIDITY.v3*GRAZ + grouping*arid_sq + arid_sq*GRAZ + (1|SITE_ID), data = twosp_dist)
Anova(twosp2)
summary(twosp2)

twosp2_simres <- simulateResiduals(twosp2)
plot(twosp2_simres)

twosp_null <- glmmTMB(euclidean_dist ~ 1 + (1|SITE_ID), data = twosp_dist)
anova(twosp_null, twosp2) #whole model p = < 2.2e-16

#get R squared
r.squaredGLMM(twosp2) #take the theoretical



#PLOTS#
#dist ~ graz
ggplot(twosp_dist, aes(x = GRAZ, y = euclidean_dist)) +
  geom_boxplot()

#dist~graz, colour by grouping
ggplot(twosp_dist, aes(x = GRAZ, y = euclidean_dist, fill = grouping)) +
  geom_boxplot()

#dist~grouping, colour by graz
ggplot(twosp_dist, aes(x = grouping, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot()


#dist~aridity
ggplot(twosp_dist , aes(x = ARIDITY.v3, y = euclidean_dist, color = grouping)) +
  geom_jitter(width = 0.1, height = 0.1)
