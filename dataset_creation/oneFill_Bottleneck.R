###############
#### SETUP ####
###############

#Load required packages
# library(Rmisc)
# library(dplyr)
# library(data.table)
# library(broom)

library(Rmisc, lib = "~/Rpackages")
library(dplyr)
library(data.table, lib = "~/Rpackages")
library(broom)


#Working directory for datasets varying migration rate and bottleneck proportion
# setwd("/Users/jamessantangelo/Desktop/CSV/raw-data/oneFill_Bottlenecks")
setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/raw-data/oneFill_Bottlenecks')

# # Globals
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
args <- list('None', 'Low', 'High')
# args <- commandArgs(trailingOnly = TRUE)
merge_lm <- list()
merge_FirstGen <- list()
num_patches <- 40

for (i in 1:length(args)){

  print(i)
  # Load dataset that varies the bottleneck proportion and add distance column
  colsToKeep <- c("x", "y","bot", "Sim", "Generation", "Cyan", "Mat_full", "Pop_size", "Mig_rate", "pA", "pB")
  colClasses <- list(numeric = c("Cyan", "bot", "Mig_rate"),
                   integer = c("x", "y", "Sim", "Generation", "Mat_full", "Pop_size"))
  name <-  sprintf('20170924_Drift_Mig-%s_Merged.csv', args[i])
  dat <- fread(name, select = colsToKeep, colClasses = colClasses, header = T)
  # print(str(dat))
  dat$Distance  <- num_patches - sqrt((dat$x - 0)^2 + (dat$y - 0)^2)
  # print(head(dat))
  dat$Mig_rate <- as.factor(as.character(dat$Mig_rate))
  dat$bot <- as.factor(as.character(dat$bot))

  #Run model testing for change in HCN frequency with distance across matrix. Performed separately for every simulation and generation, begining with the generation the matrix fill.
  dat_lm <- dat %>%
    filter(Mat_full == 1) %>%
    group_by(Sim, Mig_rate, bot, Generation) %>%
    do(FitSim = lm(Cyan ~ Distance, data = .))

  dat_FreqFirstGen <- dat %>%
    group_by(Sim, Mig_rate, bot, Distance) %>%
    slice(which.min(Generation))

  #Create data frame with results from linear morm(dels
  FitSimCoef = tidy(dat_lm, FitSim)

  #Remove initial datasets
  rm(dat, dat_lm)

  #Subset data frame to include only slopes and P-values for the effect of distance
  dataset_lm <- FitSimCoef %>%
    filter(term == "Distance") %>%
    select(estimate, p.value)

  merge_lm[[i]] <- dataset_lm
  merge_FirstGen[[i]] <- dat_FreqFirstGen

  #Write dataset with all models to csv
  # name = sprintf("FitSimCoef_Mig-%s.csv", args[i])
  # fwrite(FitSimCoef, file = paste(today, name, sep = "_"), sep = ",", col.names = TRUE)

}

# setwd("/Users/jamessantangelo/Desktop/CSV/summary-datasets/oneFill_Bottlenecks")
setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/summary-datasets/oneFill_Bottlenecks')

merged_lm <- Reduce(function(...) merge(..., all = T), merge_lm)
fwrite(merged_lm, file = paste(today, "RegSummary_oneFill_Bottlenecks.csv", sep = "_"), sep = ",", col.names = TRUE)

merged_FreqFirstGen <- Reduce(function(...) merge(..., all = T), merge_FirstGen)
fwrite(merged_FreqFirstGen, file = paste(today, "FreqFirstGen_oneFill_Bottlenecks.csv", sep = "_"), sep = ",", col.names = TRUE)

SlopeSum_Gen <- merged_lm %>%
  group_by(Sim, bot, Mig_rate) %>%
  # filter(Generation %in% seq(from = min(Generation), to = max(Generation), by = 7)) %>%
  mutate(seq = 1:n()) %>%
  group_by(bot, Mig_rate, seq) %>%
  summarize(mean = mean(estimate),
            sd = sd(estimate),
            n = length(estimate),
            se = (sd/sqrt(n)),
            ci_mean = 1.96*se,

            prop_sigPos = (sum(estimate > 0 & p.value < 0.05)/n),
            prop_pos = (sum(estimate > 0)/n),
            se_pos = sqrt((prop_pos*(1 - prop_pos)/n)),
            ci_pos = 1.96*se_pos,
            se_sigPos = sqrt((prop_sigPos*(1 - prop_sigPos)/n)),
            ci_sigPos = 1.96*se_sigPos,

            prop_sigNeg = (sum(estimate < 0 & p.value < 0.05)/n),
            prop_neg =(sum(estimate < 0)/n),
            se_neg = sqrt((prop_neg*(1 - prop_neg)/n)),
            ci_neg = 1.96*se_neg,
            se_sigNeg = sqrt((prop_sigNeg*(1 - prop_sigNeg)/n)),
            ci_sigNeg = 1.96*se_sigNeg)


fwrite(SlopeSum_Gen, file = paste(today, "MeansProps_oneFill_Bottlenecks.csv", sep = "_"), sep = ",", col.names = TRUE)

