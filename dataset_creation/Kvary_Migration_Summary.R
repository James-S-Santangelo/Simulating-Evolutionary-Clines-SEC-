#Load required packages
library(broom)
library(data.table, lib="~/Rpackages")
library(Rmisc, lib = "~/Rpackages")
library(dplyr)

setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/Drift.Migration/1D/Kvary_AllFill')

#Load dataset varying migration rate and add distance column
colsToKeep <- c("x", "y","Mig_rate", "Sim", "Generation", "pA", "pB", "Pop_size", "K", "Cyan")
dat_Kvary_MigVary <- fread('20170810_Kvary_Migration_Merged.csv', select = colsToKeep, header = T)
dat_Kvary_MigVary$Distance  <- sqrt((dat_Kvary_MigVary$x - 0)^2 + (dat_Kvary_MigVary$y - 0)^2)

#Run model testing for change in HCN, pA, and pB frequency with distance across matrix.
#Performed separately for every simulation.
dat_Kvary_Mig_lm <- dat_Kvary_MigVary %>%
  # filter(Generation == max(Generation)) %>%
  group_by(Mig_rate, Sim, Generation) %>%
  do(FitMigSimCyan = lm(Cyan ~ Distance, data = .),
     FitMigSimpA = lm(pA ~ Distance, data = .),
     FitMigSimpB = lm(pB ~ Distance, data = .))

#Create data frame with results from linear models
FitKvary_Mig_Coef_Cyan = tidy(dat_Kvary_Mig_lm, FitMigSimCyan)
FitKvary_Mig_Coef_pA = tidy(dat_Kvary_Mig_lm, FitMigSimpA)
FitKvary_Mig_Coef_pB = tidy(dat_Kvary_Mig_lm, FitMigSimpB)

#Remove initial datasets
rm(dat_Kvary_MigVary, dat_Kvary_Mig_lm)

#Function to filter summarize linear model datasets
#Subset data frame to include only slopes and P-values for the effect of distance
Sumlm <- function(df) {
  df %>%
    filter(term == "Distance") %>%
    select(estimate, p.value)
}

#Apply summary function to dataframes with estimates from models
FitKvary_Mig_Coef_Cyan <- Sumlm(FitKvary_Mig_Coef_Cyan)
FitKvary_Mig_Coef_pA <- Sumlm(FitKvary_Mig_Coef_pA)
FitKvary_Mig_Coef_pB <- Sumlm(FitKvary_Mig_Coef_pB)

#rbind above datasets with and include dataset ID
FitKvary_Mig_Coef <- bind_rows("Cyan" = FitKvary_Mig_Coef_Cyan,
                               "pA" = FitKvary_Mig_Coef_pA,
                               "pB" = FitKvary_Mig_Coef_pB,
                               .id = "id")

#Write dataset with all models to csv
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
fwrite(FitKvary_Mig_Coef, file = paste(today, "FitKvary_Mig_Coef.csv", sep = "_"), sep = ",", col.names = TRUE)

#Get mean slope and proportion of significantly positive and negative slopes
#from linear models. Averaged across simulations.
#Confidence intervals are also calculated.
Kvary_Mig_Summary <- FitKvary_Mig_Coef %>%
  group_by(Mig_rate, id) %>%
  summarise(n = length(estimate),
            mean = mean(estimate),
            sd = sd(estimate),
            se = (sd/sqrt(n)),
            ci = 1.96*se,
            prop_sigPos = (sum(estimate > 0 & p.value < 0.05)/length(estimate)),
            se_Pos = sqrt((prop_sigPos*(1 - prop_sigPos)/length(estimate))),
            ci.Pos = 1.96*se_Pos,
            prop_sigNeg = (sum(estimate < 0 & p.value < 0.05)/length(estimate)),
            se_Neg = sqrt((prop_sigNeg*(1 - prop_sigNeg)/length(estimate))),
            ci.Neg = 1.96*se_Neg)

#Wrtie dataset to disk
fwrite(Kvary_Mig_Summary, file = paste(today, "Kvary_Mig_Summary.csv", sep = "_"), sep = ",", col.names = TRUE)
