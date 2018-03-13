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
# setwd("/Users/jamessantangelo/Desktop/CSV/raw-data/oneFill_Bottlenecks")
setwd('/Users/jamessantangelo/Desktop/CSV/raw-data/oneFill_Bottlenecks')

# # Globals
today <- gsub("-","",format(Sys.Date(), formate = "$Y$m$d"))
args1 <- list('0.00', '0.01', '0.05')
args2 <- list('0.010', '0.020', '0.035', '0.050', '0.075', '0.100', '0.200', '0.500', '0.750', '1.000')
merge_lm <- list()
merge_FirstGen <- list()
merge_PopSize <- list()
num_patches <- 15

for (i in 1:length(args1)){

  for (j in 1:length(args2)){

    print(i)
    # Load dataset that varies the bottleneck proportion and add distance column
    colsToKeep <- c("x", "y","bot", "Sim", "Generation", "Cyan", "Mat_full", "Pop_size", "Mig_rate", "pA", "pB")
    colClasses <- list(numeric = c("Cyan", "bot", "Mig_rate"),
                     integer = c("x", "y", "Sim", "Generation", "Mat_full", "Pop_size"))
    name <- sprintf('oneFill_Bottlenecks(m%s)(bot%s).csv', args1[i], args2[j])
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

    #Get frequency of HCN in first generation of population's existence for every population
    dat_FreqFirstGen <- dat %>%
      group_by(Sim, Mig_rate, bot, Distance) %>%
      slice(which.min(Generation))

    #Get size of population for every population every generation, averaged across simulations
    dat_PopSize <- dat %>%
      filter(Mat_full == 1 & Distance == min(Distance)) %>%
      group_by(Sim, Mig_rate, bot, Distance) %>%
      mutate(seq = 1:n()) %>%
      select(bot, Sim, Generation, Pop_size, Distance, seq)

    #Create data frame with results from linear models
    FitSimCoef = tidy(dat_lm, FitSim)

    #Remove initial datasets
    rm(dat, dat_lm)

    #Subset data frame to include only slopes and P-values for the effect of distance
    dataset_lm <- FitSimCoef %>%
      filter(term == "Distance") %>%
      select(estimate, p.value)

    dataset <- paste("Dataset", i, j, sep="")
    merge_lm[[dataset]] <- dataset_lm
    merge_FirstGen[[dataset]] <- dat_FreqFirstGen
    merge_PopSize[[dataset]] <- dat_PopSize
  }
}

# setwd("/Users/jamessantangelo/Desktop/CSV/summary-datasets/oneFill_Bottlenecks")
setwd('/Users/jamessantangelo/Desktop/CSV/summary-datasets/oneFill_Bottlenecks')

merged_lm <- Reduce(function(...) merge(..., all = T), merge_lm)
fwrite(merged_lm, file = paste(today, "RegSummary_oneFill_Bottlenecks.csv", sep = "_"), sep = ",", col.names = TRUE)

merged_FreqFirstGen <- Reduce(function(...) merge(..., all = T), merge_FirstGen)
fwrite(merged_FreqFirstGen, file = paste(today, "FreqFirstGen_oneFill_Bottlenecks.csv", sep = "_"), sep = ",", col.names = TRUE)

merged_PopSize <- Reduce(function(...) merge(..., all = T), merge_PopSize)
fwrite(merged_PopSize, file = paste(today, "PopSize_oneFill_Bottlenecks.csv", sep = "_"), sep = ",", col.names = TRUE)

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

