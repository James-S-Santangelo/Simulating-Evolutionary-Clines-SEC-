#Load required packages
library(broom)
library(data.table)
library(Rmisc)
library(dplyr)

#Load required packages
# library(broom)
# library(data.table, lib="~/Rpackages")
# library(Rmisc, lib = "~/Rpackages")
# library(dplyr)

# setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/drift-migration/1D/Kvary_AllFill')

setwd('/Users/jamessantangelo/Desktop/CSV/Kvary')

# Globals
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
args <- list('0.0000', '0.0010', '0.0025',
             '0.0050', '0.0100', '0.0200',
             '0.0350', '0.0500', '0.1000',
             '0.2000', '0.3500', '0.5000'
             '1.0000')
merge_lm <- list()

for(i in 1:length(args)){

  print(i, args[i])

  #Load dataset varying migration rate and add distance column
  colsToKeep <- c("x", "y","Mig_rate", "Sim", "Generation", "pA", "pB", "Pop_size", "K", "Cyan")
  colClasses <- list(numeric = c("Cyan", "pA", "pB", "Mig_rate"),
                    integer = c("x", "y", "Sim", "Generation", "Pop_size", "K"))
  name <-  sprintf('20171025_drift_Kvary-AllFill(m%s).csv', args[i])
  dat_Kvary_MigVary <- fread(name, select = colsToKeep, colClasses = colClasses, header = T)
  dat_Kvary_MigVary$Distance  <- sqrt((dat_Kvary_MigVary$x - 0)^2 + (dat_Kvary_MigVary$y - 0)^2)
  dat_Kvary_MigVary$Mig_rate <- as.factor(as.character(dat_Kvary_MigVary$Mig_rate))


  #Run model testing for change in HCN frequency with distance across matrix. Performed separately for every simulation and generation.
  dat_Kvary_Mig_lm <- dat_Kvary_MigVary %>%
    # filter(Generation == max(Generation)) %>%
    group_by(Mig_rate, Sim, Generation) %>%
    do(FitMigSimCyan = lm(Cyan ~ Distance, data = .),
       FitMigSimpA = lm(pA ~ Distance, data = .),
       FitMigSimpB = lm(pB ~ Distance, data = .))

  #Create data frame with results from linear models of Cyan
  FitKvary_Mig_Coef_Cyan = tidy(dat_Kvary_Mig_lm, FitMigSimCyan)
  FitKvary_Mig_Coef_pA = tidy(dat_Kvary_Mig_lm, FitMigSimpA)
  FitKvary_Mig_Coef_pB = tidy(dat_Kvary_Mig_lm, FitMigSimpB)

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

  rm(FitKvary_Mig_Coef_Cyan, FitKvary_Mig_Coef_pA, FitKvary_Mig_Coef_pB)

  merge_lm[[i]] <- FitKvary_Mig_Coef
}

#Write dataset with all models to csv
merged_lm <- Reduce(function(...) merge(..., all = T), merge_lm)
fwrite(merged_lm, file = paste(today, "Kvary-WithMig_Coef.csv", sep = "_"), sep = ",", col.names = TRUE)


#Get mean slope and proportion of significantly positive and negative slopes
#from linear models. Averaged across simulations.
#Confidence intervals are also calculated.
SlopeSum_Gen <- merged_lm %>%
  group_by(id, Mig_rate) %>%
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

#Wrtie dataset to disk
fwrite(SlopeSum_Gen, file = paste(today, "Kvary-WithMig_SlopeSummary.csv", sep = "_"), sep = ",", col.names = TRUE)
