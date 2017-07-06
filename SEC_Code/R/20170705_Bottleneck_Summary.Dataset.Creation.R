#Load required packages
library(dplyr)
library(broom)
library("data.table", lib="~/Rpackages")

#Working directory for datasets varying migration rate and bottleneck proportion
setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/Drift.Migration/1D/Mig_Bot_Vary')

# Load dataset that varies the bottleneck proportion and add distance column
colsToKeep <- c("x", "y","bot", "Sim", "Generation", "Acyan", "Mat.full")
dat_Bot_Vary <- fread('20170704_Merged_BotOnly.csv', select = colsToKeep, header = T)
dat_Bot_Vary$Distance  <- sqrt((dat_Bot_Vary$x - 0)^2 + (dat_Bot_Vary$y - 0)^2)
dat_Bot_Vary$Cyan  <- 1 - dat_Bot_Vary$Acyan

#Run model testing for change in HCN frequency with distance across matrix. 
#Performed separately for every simulation and generation, begining with the generation the matrix fill.
dat_Botlm_sum <- dat_Bot_Vary %>%
  group_by(bot, Sim, Generation) %>% 
  filter(Mat.full == 1) %>%
  do(FitBotSim = lm(Cyan ~ Distance, data = .))

#Create data frame with results from linear models
FitBotSimCoef = tidy(dat_Botlm_sum, FitBotSim)

#Remove initial datasets
rm(dat_Bot_Vary, dat_Botlm_sum)

#Subset data frame to include only slopes and P-values for the effect of distance
FitBotSimCoef <- FitBotSimCoef %>% 
  filter(term == "Distance") %>%
  select(estimate, p.value)

#Write dataset with all models to csv
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
fwrite(FitBotSimCoef, file = paste(today, "FitBotSimCoef.csv", sep = "_"), sep = ",", col.names = TRUE)

#Get mean slope and proportion of significantly positive and negative slopes
#from linear models. Done for each generation, averaged across simulations. 
#Confidence intervals are also calculated.
Bot_Vary_SlopeSum <- FitBotSimCoef %>%
  group_by(bot, Generation) %>%
  summarise(mean = mean(estimate), 
            sd = sd(estimate), 
            n = length(estimate), 
            se = (sd/sqrt(n)), 
            ci.lower = mean - 1.96*se,
            ci.upper = mean + 1.96*se,
            prop_sigPos = (sum(estimate > 0 & p.value < 0.05)/length(estimate)), 
            se_Pos = sqrt((prop_sigPos*(1 - prop_sigPos)/length(estimate))),
            ci.lower.Pos = prop_sigPos - 1.96*se_Pos,
            ci.upper.Pos = prop_sigPos + 1.96*se_Pos,
            prop_sigNeg = (sum(estimate < 0 & p.value < 0.05)/length(estimate)),
            se_Neg = sqrt((prop_sigNeg*(1 - prop_sigNeg)/length(estimate))),
            ci.lower.Neg = prop_sigNeg - 1.96*se_Neg,
            ci.upper.Neg = prop_sigNeg + 1.96*se_Neg)

#Write dataset with summary info to csv
fwrite(FitBotSimCoef, file = paste(today, "Bot_Vary_SlopeSum.csv", sep = "_"), sep = ",", col.names = TRUE)
