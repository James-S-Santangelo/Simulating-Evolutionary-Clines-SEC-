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

setwd('/Users/jamessantangelo/Desktop/CSV/raw-data/allFill_Kvary/allFill_Kvary_KminVary')
# setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/raw-data/allFill_Kvary/allFill_Kvary_KminVary')

# Globals
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
# args <- list("0.00", "0.01", "0.05")
args1 <- list('0.00', '0.01', '0.05')
args2 <- list('10', '100', '500', '1000')
merge_lm <- list()
num_patches <- 15

for(i in 1:length(args1)){

  for (j in 1:length(args2)){

    print(c(i, args1[i], j, args2[j]))

    #Load dataset varying migration rate and add distance column
    colsToKeep <- c("x", "y","Mig_rate", "Sim", "Generation", "pA", "pB", "Pop_size", "K", "Cyan", "min_K")
    colClasses <- list(numeric = c("Cyan", "pA", "pB", "Mig_rate"),
                      integer = c("x", "y", "Sim", "Generation", "Pop_size", "K", "min_K"))
    name <-  sprintf('allFill_Kvary(Kmin%s)(m%s).csv', args2[j], args1[i])
    dat_Kvary_MigVary <- fread(name, select = colsToKeep, colClasses = colClasses, header = T)
    dat_Kvary_MigVary$Distance  <- num_patches - sqrt((dat_Kvary_MigVary$x - 0)^2 + (dat_Kvary_MigVary$y - 0)^2)
    dat_Kvary_MigVary$Mig_rate <- as.factor(as.character(dat_Kvary_MigVary$Mig_rate))
    dat_Kvary_MigVary$min_K <- as.factor(as.character(dat_Kvary_MigVary$min_K))


    #Run model testing for change in HCN frequency with distance across matrix. Performed separately for every simulation and generation.
    range <- max(dat_Kvary_MigVary$Distance) - min(dat_Kvary_MigVary$Distance)
    dat_Kvary_Mig_lm <- dat_Kvary_MigVary %>%
      # filter(Generation == max(Generation)) %>%
      group_by(Mig_rate, min_K, Sim, Generation) %>%
      mutate(Distance_std = ((Distance - min(Distance)) / range)) %>%
      do(FitMigSimCyan = lm(Cyan ~ Distance_std, data = .))

    #Create data frame with results from linear models of Cyan
    FitKvary_Mig_Coef_Cyan = tidy(dat_Kvary_Mig_lm, FitMigSimCyan)

    rm(dat_Kvary_MigVary, dat_Kvary_Mig_lm)

    #Function to filter summarize linear model datasets
    #Subset data frame to include only slopes and P-values for the effect of distance
    dataset_lm <- FitKvary_Mig_Coef_Cyan %>%
      filter(term == "Distance_std") %>%
      select(estimate, p.value)

    dataset <- paste("Dataset", i, j, sep="")
    merge_lm[[dataset]] <- dataset_lm
  }
}

setwd('/Users/jamessantangelo/Desktop/CSV/summary-datasets/allFill_Kvary/allFill_Kvary_KminVary')
# setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/summary-datasets/allFill_Kvary/allFill_Kvary_KminVary')

#Write dataset with all models to csv
merged_lm <- Reduce(function(...) merge(..., all = T), merge_lm)
fwrite(merged_lm, file = paste(today, "StdSlopes_allFill_Kvary_KminVary.csv", sep = "_"), sep = ",", col.names = TRUE)

