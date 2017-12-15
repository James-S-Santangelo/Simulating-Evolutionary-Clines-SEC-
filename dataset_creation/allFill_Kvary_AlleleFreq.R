###############
#### SETUP ####
###############

#Load required packages
library(Rmisc)
library(dplyr)
library(data.table)
library(broom)

# library(Rmisc, lib = "~/Rpackages")
# library(dplyr)
# library(data.table, lib = "~/Rpackages")
# library(broom)


#Working directory for datasets varying migration rate and bottleneck proportion
setwd("/Users/jamessantangelo/Desktop/CSV/allFill_Kvary_AlleleFreq")
# setwd('/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/Drift.Migration/1D/Mig_Bot_Vary')


# # Globals
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
args1 <- list('0.10', '0.50', '0.90')
args2 <- list('0.10', '0.50', '0.90')
args3 <- list('0.00', '0.01', '0.05')

merge_lm = list()
num_patches <- 40

# merge_lm <- list()
# merge_FirstGen <- list()

for (i in 1:length(args1)){

  for (j in 1:length(args2)){

    for (k in 1:length(args3)){


      print(c(i, j, k))
      # Load dataset that varies the bottleneck proportion and add distance column
      colsToKeep <- c("x", "y","bot", "Sim", "Generation", "Cyan", "Mat_full", "Pop_size", "Mig_rate", "pA", "pB")
      colClasses <- list(numeric = c("Cyan", "bot", "Mig_rate"),
                       integer = c("x", "y", "Sim", "Generation", "Mat_full", "Pop_size"))
      name <-  sprintf('allFill_Kvary_AlleleFreq(pA%s)(pB%s)(m%s).csv', args1[i], args2[j], args3[k])
      dat <- fread(name, select = colsToKeep, colClasses = colClasses, header = T)
      # print(str(dat))
      dat$Distance  <- num_patches - sqrt((dat$x - 0)^2 + (dat$y - 0)^2)
      # print(head(dat))
      dat$pA <- as.factor(as.character(dat$pA))
      dat$pB <- as.factor(as.character(dat$pB))
      dat$Mig_rate <- as.factor(as.character(dat$Mig_rate))

      #Add starting pA and pB as column to dataframe
      dat <- dat %>%
        mutate(pA_start = args1[[i]], pB_start = args2[[j]])

      #Run model testing for change in HCN frequency with distance across matrix. Performed separately for every simulation and generation, begining with the generation the matrix fill.
      dat_lm <- dat %>%
        filter(Mat_full == 1) %>%
        group_by(Sim, Mig_rate, pA_start, pB_start, Generation) %>%
        do(FitSim = lm(Cyan ~ Distance, data = .))

      #Create data frame with results from linear morm(dels
      FitSimCoef = tidy(dat_lm, FitSim)

      #Remove initial datasets
      rm(dat, dat_lm)

      #Subset data frame to include only slopes and P-values for the effect of distance
      dataset_lm <- FitSimCoef %>%
        filter(term == "Distance") %>%
        select(estimate, p.value)

      # print(head(dataset_lm))
      dataset <- paste("Dataset", i, j, k, sep="")
      merge_lm[[dataset]] <- dataset_lm
    }
  }
}

merged_lm <- Reduce(function(...) merge(..., all = T), merge_lm)
fwrite(merged_lm, file = paste(today, "RegSummary_allFill_Kvary_AlleleFreq.csv", sep = "_"), sep = ",", col.names = TRUE)

SlopeSum_Gen <- merged_lm %>%
  group_by(Sim, Mig_rate, pA_start, pB_start) %>%
  # filter(Generation %in% seq(from = min(Generation), to = max(Generation), by = 7)) %>%
  mutate(seq = 1:n()) %>%
  group_by(Mig_rate, pA_start, pB_start, seq) %>%
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


fwrite(SlopeSum_Gen, file = paste(today, "MeansProps_allFill_Kvary_AlleleFreq.csv", sep = "_"), sep = ",", col.names = TRUE)






