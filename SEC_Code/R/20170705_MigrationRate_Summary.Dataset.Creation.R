#Load required packages
library(broom)
library("data.table", lib="~/Rpackages")
library(Rmisc, lib = "~/Rpackages")
library(dplyr)

#Working directory for datasets varying migration rate and bottleneck proportion
setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/Drift.Migration/1D/Mig_Bot_Vary')

#Load dataset varying migration rate and add distance column
colsToKeep <- c("x", "y", "Mig_rate", "Sim", "Generation", "Cyan", "Mat_full", "Pop_size")
dat_Mig_Vary <- fread('20170912_Merged_MigOnly.csv', select = colsToKeep, header = T)
dat_Mig_Vary$Distance  <- sqrt((dat_Mig_Vary$x - 0)^2 + (dat_Mig_Vary$y - 0)^2)
# dat_Mig_Vary$Recip = 1 / dat_Mig_Vary$Pop_size

#Generate dataset showing Ne for every population, grouped by bottleneck strength
# dat_mig_Ne <- dat_Mig_Vary %>%
#   group_by(Sim, Population, Mig_rate, Distance) %>%
#   summarise(sumRecip = sum(Recip),
#             Generations = length(unique(Generation)),
#             Ne = Generations / sumRecip)
# dat_mig_Ne <- summarySE(dat_mig_Ne, groupvars = c("Population", "Mig_rate", "Distance"), measurevar = "Ne")

#Write Ne dataset to csv
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
# fwrite(dat_mig_Ne, file = paste(today, "Ne_Mig.csv", sep = "_"), sep = ",", col.names = TRUE)

#Run model testing for change in HCN frequency with distance across matrix.
#Performed separately for every simulation and generation, begining with the generation the matrix full.
dat_Miglm_sum <- dat_Mig_Vary %>%
  group_by(Mig_rate, Sim, Generation) %>%
  filter(Mat_full == 1) %>%
  do(FitMigSim = lm(Cyan ~ Distance, data = .))

#Create data frame with results from linear models
FitMigSimCoef = tidy(dat_Miglm_sum, FitMigSim)

#Remove initial datasets
rm(dat_Mig_Vary, dat_Miglm_sum)

#Subset data frame to include only slopes and P-values for the effect of distance
FitMigSimCoef <- FitMigSimCoef %>%
  filter(term == "Distance") %>%
  select(estimate, p.value)

#Write dataset with all models to csv
fwrite(FitMigSimCoef, file = paste(today, "FitMigSimCoef.csv", sep = "_"), sep = ",", col.names = TRUE)

#Get mean slope and proportion of significantly positive and negative slopes
#from linear models. Done for each generation, averaged across simulations.
#Confidence intervals are also calculated.
MigRate_SlopeSum_Gen <- FitMigSimCoef %>%
  group_by(Mig_rate, Generation) %>%
  summarise(mean = mean(estimate),
            sd = sd(estimate),
            n = length(estimate),
            se = (sd/sqrt(n)),
            ci.lower = 1.96*se,
            ci.upper = 1.96*se,
            prop_sigPos = (sum(estimate > 0 & p.value < 0.05)/length(estimate)),
            se_Pos = sqrt((prop_sigPos*(1 - prop_sigPos)/length(estimate))),
            ci.lower.Pos = 1.96*se_Pos,
            ci.upper.Pos = 1.96*se_Pos,
            prop_sigNeg = (sum(estimate < 0 & p.value < 0.05)/length(estimate)),
            se_Neg = sqrt((prop_sigNeg*(1 - prop_sigNeg)/length(estimate))),
            ci.lower.Neg = 1.96*se_Neg,
            ci.upper.Neg = 1.96*se_Neg)

#Write dataset with summary info to csv
fwrite(MigRate_SlopeSum_Gen, file = paste(today, "MigRate_SlopeSum_Gen.csv", sep = "_"), sep = ",", col.names = TRUE)

#Get mean slope and proportion of significantly positive and negative slopes
#for each value of the migration rate. 95% CI's also calculated.
MigRate_SlopeSum <- FitMigSimCoef %>%
  group_by(Mig_rate) %>%
  summarise(mean = mean(estimate),
            sd = sd(estimate),
            n = length(estimate),
            se = (sd/sqrt(n)),
            ci.lower = 1.96*se,
            ci.upper = 1.96*se,
            prop_sigPos = (sum(estimate > 0 & p.value < 0.05)/length(estimate)),
            se_Pos = sqrt((prop_sigPos*(1 - prop_sigPos)/length(estimate))),
            ci.lower.Pos = 1.96*se_Pos,
            ci.upper.Pos = 1.96*se_Pos,
            prop_sigNeg = (sum(estimate < 0 & p.value < 0.05)/length(estimate)),
            se_Neg = sqrt((prop_sigNeg*(1 - prop_sigNeg)/length(estimate))),
            ci.lower.Neg = 1.96*se_Neg,
            ci.upper.Neg = 1.96*se_Neg)

#Wrtie dataset to disk
fwrite(MigRate_SlopeSum, file = paste(today, "MigRate_SlopeSum.csv", sep = "_"), sep = ",", col.names = TRUE)
