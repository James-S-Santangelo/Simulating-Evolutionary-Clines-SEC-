###############
#### SETUP ####
###############

#Load required packages
library(Rmisc)
library(dplyr)
# library(dplyr, lib = "~/Rpackages")
library(data.table)
# library(data.table, lib = "~/Rpackages")
library(broom)

#Working directory for datasets varying migration rate and bottleneck proportion
setwd("/Users/jamessantangelo/Desktop/CSV")
# setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/Drift.Migration/1D/Mig_Bot_Vary')

# # Globals
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
args <- list('None', 'Low', 'High')
# args <- commandArgs(trailingOnly = TRUE)
datasets_to_merge <- list()

for (i in 1:length(args)){

  print(i)
  # Load dataset that varies the bottleneck proportion and add distance column
  colsToKeep <- c("x", "y","bot", "Sim", "Generation", "Cyan", "Mat_full", "Pop_size", "Mig_rate")
  colClasses <- list(numeric = c("Cyan"),
                   factor = c("bot", "Mig_rate"),
                   integer = c("x", "y", "Sim", "Generation", "Mat_full", "Pop_size"))
  name <-  sprintf('20170923_Drift_Mig-%s_Merged_R-test.csv', args[i])
  dat <- fread(name, select = colsToKeep, colClasses = colClasses, header = T)
  # print(str(dat))
  dat$Distance  <- sqrt((dat$x - 0)^2 + (dat$y - 0)^2)
  # print(head(dat))

  #Run model testing for change in HCN frequency with distance across matrix. Performed separately for every simulation and generation, begining with the generation the matrix fill.
  dat_lm <- dat %>%
    filter(Mat_full == 1) %>%
    group_by(Sim, Mig_rate, bot, Generation) %>%
    do(FitSim = lm(Cyan ~ Distance, data = .))

  #Create data frame with results from linear morm(dels
  FitSimCoef = tidy(dat_lm, FitSim)

  #Remove initial datasets
  rm(dat, dat_lm)

  #Subset data frame to include only slopes and P-values for the effect of distance
  dataset = sprintf("FitSimCoef_Mig-%s.csv", args[i])
  dataset <- FitSimCoef %>%
    filter(term == "Distance") %>%
    select(estimate, p.value)

  datasets_to_merge[[i]] <- dataset

  #Write dataset with all models to csv
  # name = sprintf("FitSimCoef_Mig-%s.csv", args[i])
  # fwrite(FitSimCoef, file = paste(today, name, sep = "_"), sep = ",", col.names = TRUE)

}

merged <- Reduce(function(...) merge(..., all = T), datasets_to_merge)
# fwrite(merged, file = paste(today, "FitSimCoef_Mig-Merged.csv", sep = "_"), sep = ",", col.names = TRUE)

SlopeSum_Gen <- merged %>%
  group_by(Sim, bot, Mig_rate) %>%
  # filter(Generation %in% seq(from = min(Generation), to = max(Generation), by = 7)) %>%
  mutate(seq = 1:n()) %>%
  group_by(bot, Mig_rate, seq) %>%
  summarize(mean = mean(estimate),
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
  




# mutate(group_by(filtered, Sim, bot, Mig_rate), seq = 1:n())
# summarise(mean = mean(estimate))

            # sd = sd(estimate),
            # n = length(estimate),
            # se = (sd/sqrt(n)),
            # ci.lower = 1.96*se,
            # ci.upper = 1.96*se,
            # prop_sigPos = (sum(estimate > 0 & p.value < 0.05)/length(estimate)),
            # se_Pos = sqrt((prop_sigPos*(1 - prop_sigPos)/length(estimate))),
            # ci.lower.Pos = 1.96*se_Pos,
            # ci.upper.Pos = 1.96*se_Pos,
            # prop_sigNeg = (sum(estimate < 0 & p.value < 0.05)/length(estimate)),
            # se_Neg = sqrt((prop_sigNeg*(1 - prop_sigNeg)/length(estimate))),
            # ci.lower.Neg = 1.96*se_Neg,
            # ci.upper.Neg = 1.96*se_Neg)

# test

  # group_by(bot, Mig_rate) #%>%
  # summarize(mean = mean(estimate))

#Get mean slope and proportion of significantly positive and negative slopes from linear models. Done for each generation, averaged across simulations. Confidence intervals are also calculated.
# SlopeSum_Gen <- merged %>%
#   group_by(bot, Mig_rate, Generation) %>%
#   summarise(mean = mean(estimate),
#             sd = sd(estimate),
#             n = length(estimate),
#             se = (sd/sqrt(n)),
#             ci.lower = 1.96*se,
#             ci.upper = 1.96*se,
#             prop_sigPos = (sum(estimate > 0 & p.value < 0.05)/length(estimate)),
#             se_Pos = sqrt((prop_sigPos*(1 - prop_sigPos)/length(estimate))),
#             ci.lower.Pos = 1.96*se_Pos,
#             ci.upper.Pos = 1.96*se_Pos,
#             prop_sigNeg = (sum(estimate < 0 & p.value < 0.05)/length(estimate)),
#             se_Neg = sqrt((prop_sigNeg*(1 - prop_sigNeg)/length(estimate))),
#             ci.lower.Neg = 1.96*se_Neg,
#             ci.upper.Neg = 1.96*se_Neg)

# #Write dataset with summary info to csv
# fwrite(Bot_SlopeSum_Gen, file = paste(today, "Bot_SlopeSum_Gen.csv", sep = "_"), sep = ",", col.names = TRUE)

# #Get mean slope and proportion of significantly positive and negative slopes
# #for each value of the bottleneck proportion. 95% CI's also calculated.
# Bot_SlopeSum <- FitBotSimCoef %>%
#   group_by(bot) %>%
#   summarise(mean = mean(estimate),
#             sd = sd(estimate),
#             n = length(estimate),
#             se = (sd/sqrt(n)),
#             ci.lower = 1.96*se,
#             ci.upper = 1.96*se,
#             prop_sigPos = (sum(estimate > 0 & p.value < 0.05)/length(estimate)),
#             se_Pos = sqrt((prop_sigPos*(1 - prop_sigPos)/length(estimate))),
#             ci.lower.Pos = 1.96*se_Pos,
#             ci.upper.Pos = 1.96*se_Pos,
#             prop_sigNeg = (sum(estimate < 0 & p.value < 0.05)/length(estimate)),
#             se_Neg = sqrt((prop_sigNeg*(1 - prop_sigNeg)/length(estimate))),
#             ci.lower.Neg = 1.96*se_Neg,
#             ci.upper.Neg = 1.96*se_Neg)
#
# #Write dataset with summary info to csv
# fwrite(Bot_SlopeSum, file = paste(today, "Bot_SlopeSum.csv", sep = "_"), sep = ",", col.names = TRUE)
#

#### EFFECTIVE POPULATION SIZE ####

# dat_Bot_Vary$Recip = 1 / dat_Bot_Vary$Pop_size

#Generate dataset showing Ne for every population, grouped by bottleneck strength
# dat_bot_Ne <- dat_Bot_Vary %>%
#   group_by(Sim, Population, bot, Distance) %>%
#   summarise(sumRecip = sum(Recip),
#             Generations = length(unique(Generation)),
#             Ne = Generations / sumRecip)
# dat_bot_Ne <- summarySE(dat_bot_Ne, groupvars = c("Population", "bot", "Distance"), measurevar = "Ne")

#Write Ne dataset to csv
# today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
# fwrite(dat_bot_Ne, file = paste(today, "Ne_Bot.csv", sep = "_"), sep = ",", col.names = TRUE)
