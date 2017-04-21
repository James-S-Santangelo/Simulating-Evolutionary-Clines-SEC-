#Load required packages
library(data.table)
library(plyr)
library(ggplot2)

#Dataset varying population size
setwd("/Users/jamessantangelo/Documents/Academia/Master's/SEC - Simulating evolutionary clines/Drift/Data, R code and figures/Datasets/Nvary")
datNvary<-fread("20160906_SEC_Drift_Nvary.results.csv",header = T)

#Dataset varying number of generations
setwd("/Users/jamessantangelo/Documents/Academia/Master's/SEC - Simulating evolutionary clines/Drift/Data, R code and figures/Datasets/StepVary")
datStepVary<-fread("20160907_SEC_Drift_StepVary.results.csv",header = T)

#Dataset varying allele frequencies
setwd("/Users/jamessantangelo/Documents/Academia/Master's/SEC - Simulating evolutionary clines/Drift/Data, R code and figures/Datasets/pApBvary")
datpApBVary<-fread("20160906_SEC_Drift_pApBvary.results.csv",header = T)
datpApBVary$pA.pB<-paste(datpApBVary$pAi,datpApBVary$pBi,sep=";")


#Theme used to plot figures throughout script. Modify depending on usage (e.g. poster, talk, etc.)
ng1=theme(aspect.ratio=0.7,panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line.x = element_line(color="black",size=1), 
          axis.line.y = element_line(color="black",size=1),
          axis.ticks=element_line(color="black"), 
          axis.text=element_text(color="black",size=15), 
          axis.title=element_text(color="black",size=1), 
          axis.title.y=element_text(vjust=2,face="bold",size=15),
          axis.title.x=element_text(vjust=0.1,face="bold",size=15),
          axis.text.x=element_text(size=13),
          axis.text.y=element_text(size=13),
          legend.position = "right", legend.direction="vertical", 
          legend.text=element_text(size=11), legend.key = element_rect(fill = "white"), 
          legend.title = element_text(size=13,face="bold"),legend.key.size = unit(0.5, "cm"))

##############################
#### VARY POPULATION SIZE ####
##############################

##Summary dataset showing mean phenotype frequency (recessive) and mean allele frequencies
#by generation number and population size. Averaged across number of simulations
NvarySum<-ddply(datNvary,.(N,step),summarise,Phen = mean(Phen),pA = mean(pA),pB = mean(pB))

##Plot of the frequency of recessive phenotype by the number of generations with 
#different lines for different population sizes. Averaged across 1000 simulations
plotPhen.StepNvary<-ggplot(NvarySum, aes(x = step, y = Phen,group = N)) + 
  ylab("Frequency of recessive
phenotype")+xlab("Generation number")+geom_line(size=0.5, aes(color = factor(N)))+
  labs(colour = "Pop. size")+
  coord_cartesian(ylim = c(0.42,0.78))+
  scale_y_continuous(breaks=seq(from = 0.42, to = 0.78, by = 0.04))
plotPhen.StepNvary+ng1

#Dataset containing the first row from each simulation where the recessive phenotype is fixed
Nvary_Phen.fixU <- datNvary[datNvary$Phen == 1.0,]
Nvary_Phen.fixU <- Nvary_Phen.fixU[,unique(Nvary_Phen.fixU, by = c('N','sim'))]

#Dataset containing the first row from each simulation where the dominant phenotype is fixed
Nvary_Phen.fixD <- datNvary[datNvary$Phen == 0.0,]
Nvary_Phen.fixD <- Nvary_Phen.fixD[,unique(Nvary_Phen.fixD, by = c('N','sim'))]

#Datasets for number of simulations resulting in fixation of either phenotpype
Nvary_fixU <- ddply(Nvary_Phen.fixU,.(N),summarise,prop=(length(N)/1000))
Nvary_fixU <- rbind(Nvary_fixU,c(500,0.000))
Nvary_fixD <- ddply(Nvary_Phen.fixD,.(N),summarise,prop=(length(N)/1000))
Nvary_fixD <- rbind(Nvary_fixD,c(200,0.000),c(300,0.000),c(400,0.000),c(500,000))
Nvary_fix <- merge(Nvary_fixU,Nvary_fixD,by.x ="N",by.y ="N")
setnames(Nvary_fix,old=c("prop.x","prop.y"),new=c("FixU","FixD"))
Nvary_fix$FixTot <- Nvary_fix$FixU + Nvary_fix$FixD

##Dataset showing mean time to fixation of recessive phenotype for each population size
#Also mean frequency of both alleles at phenotype fixation
Nvary_fixGen <- ddply(Nvary_Phen.fixU,.(N),summarise,meanGen = mean(step),
                      mean.pA = mean(pA),mean.pB = mean(pB))

#Plot showing mean number of generation to fixation of recessive phenotype by populations size
plotGen.N.Fix<-ggplot(Nvary_fixGen, aes(x = factor(N), y = meanGen)) + 
  ylab("Mean generations
to fixation")+xlab("Population size")+geom_bar(stat="identity", color = "black",fill="grey")+
  coord_cartesian(ylim = c(0,50))+
  scale_y_continuous(breaks=seq(from = 0, to = 50, by = 5))
plotGen.N.Fix+ng1.45

####################################
#### VARY NUMBER OF GENERATIONS ####
####################################

#Theme used to plot figures throughout script. Modify depending on usage (e.g. poster, talk, etc.)
ng1.45=theme(aspect.ratio=0.7,panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line.x = element_line(color="black",size=1), 
          axis.line.y = element_line(color="black",size=1),
          axis.ticks=element_line(color="black"), 
          axis.text=element_text(color="black",size=15), 
          axis.title=element_text(color="black",size=1), 
          axis.title.y=element_text(vjust=2,face="bold",size=15),
          axis.title.x=element_text(vjust=0.1,face="bold",size=15),
          axis.text.x=element_text(size=13,angle=45,hjust=1),
          axis.text.y=element_text(size=13),
          legend.position = "right", legend.direction="vertical", 
          legend.text=element_text(size=11), legend.key = element_rect(fill = "white"), 
          legend.title = element_text(size=13,face="bold"),legend.key.size = unit(0.5, "cm"))

##Summary dataset showing mean phenotype frequency (recessive) and mean allele frequencies
#by number of generation per simulation. Averaged across 1000 simulations
StepvarySum<-ddply(datStepVary,.(Stepi),summarise,Phen = mean(Phen),pA = mean(pA),pB = mean(pB))

##Bar graph showing the mean phenotype frequency (recessive) by the number of generations
#during cline model. Averaged across all simulations
plotPhen.Stepi<-ggplot(StepvarySum, aes(x = factor(Stepi), y = Phen)) + 
  ylab("Frequency of recessive
phenotype")+xlab("Number of generations")+geom_bar(stat="identity", color = "black",fill="grey")+
  labs(colour = "Pop. size")+
  coord_cartesian(ylim = c(0.42,0.78))+
  scale_y_continuous(breaks=seq(from = 0.42, to = 0.78, by = 0.04))
plotPhen.Stepi+ng1.45

#Dataset containing the first row from each simulation where the recessive phenotype is fixed
Step_Phen.fixU <- datStepVary[datStepVary$Phen == 1.0,]
Step_Phen.fixU <- Step_Phen.fixU[,unique(Step_Phen.fixU, by = c('Stepi','sim'))]

#Dataset containing the first row from each simulation where the dominant phenotype is fixed
Step_Phen.fixD <- datStepVary[datStepVary$Phen == 0.0,]
Step_Phen.fixD <- Step_Phen.fixD[,unique(Step_Phen.fixD, by = c('Stepi','sim'))]

#Datasets for number of simulations resulting in fixation of either phenotpype
Step_fixU <- ddply(Step_Phen.fixU,.(Stepi),summarise,prop=(length(Stepi)/1000))
Step_fixU <- rbind(Step_fixU,c(2,0.000),c(4,0.000),c(6,0.000),c(8,0.000),c(10,0.000))
Step_fixD <- ddply(Step_Phen.fixD,.(Stepi),summarise,prop=(length(Stepi)/1000))
Step_fixD <- rbind(Step_fixD,c(2,0.000),c(4,0.000),c(6,0.000),c(8,0.000),c(10,0.000),c(20,0.000),c(30,0.000),c(40,0.000))
Step_fix <- merge(Step_fixU,Step_fixD,by.x ="Stepi",by.y ="Stepi")
setnames(Step_fix,old=c("prop.x","prop.y"),new=c("FixU","FixD"))
Step_fix$FixTot <- Step_fix$FixU + Step_fix$FixD

#Plot showing proportion of simulations resulting in fixation of either phenotype
plotGen.Step.Fix<-ggplot(Step_fix, aes(x = factor(Stepi), y = FixD)) + 
  ylab("Proportion fixed")+xlab("Number of generations")+geom_bar(stat="identity", color = "black",fill="grey")+
  coord_cartesian(ylim = c(0,1.0))+
  scale_y_continuous(breaks=seq(from = 0, to = 1.0, by = 0.1))
plotGen.Step.Fix+ng1.45

#########################################
#### VARY INITIAL ALLELE FREQUENCIES ####
#########################################

##Summary dataset showing mean phenotype frequency (recessive) and mean allele frequencies
#by number of generation per simulation. Averaged across number of simulations
pApBvarySum<-ddply(datpApBVary,.(pAi,pBi,step),summarise,Phen = mean(Phen),pA = mean(pA),pB = mean(pB))

###Plot of the frequency of recessive phenotype by the number of generations with 
##different lines for different starting allele frequencies. Averaged across 1000 simulations
#Not the easiest to interpret. Will produce two separate figures
plotPhen.Step.pApBvary<-ggplot(pApBvarySum, aes(x = step, y = Phen,group = pA.pB)) + 
  ylab("Frequency of recessive
phenotype")+xlab("Generation number")+geom_line(size=0.5, aes(color = factor(pA.pB)))+
  labs(colour = "pA;pB")+
  coord_cartesian(ylim=c(0,1.0))+
  scale_y_continuous(breaks=seq(from = 0, to = 1.0, by = 0.1))
plotPhen.Step.pApBvary+ng1

##Create two dataset: one for difference in initial allele frequencies and one for
#indentical inutial allele frequencies.
pApBvarySum_Diff<-pApBvarySum[pApBvarySum$pAi != pApBvarySum$pBi,]
pApBvarySum_Same<-pApBvarySum[pApBvarySum$pAi == pApBvarySum$pBi,]

#Calcuale Difference in initial allele frequencies.
pApBvarySum_Diff$Diff = pApBvarySum_Diff$pBi - pApBvarySum_Diff$pAi

##Plot frequency of recessive phenotype by generation number with different lines
#for differences in initial allele frequencies.
plotPhen.pApBdiff<-ggplot(pApBvarySum_Diff, aes(x = step, y = Phen, group = Diff)) + 
  ylab("Frequency of recessive
phenotype")+xlab("Generation number")+
  geom_line(size = 0.5,aes(color = factor(Diff)))+
  coord_cartesian(ylim = c(0.46,0.86))+
  scale_y_continuous(breaks=seq(from = 0.46, to = 0.86, by = 0.04))+
  labs(colour = "Difference 
(pB - pA)")+
  guides(colour=guide_legend(reverse=TRUE))
plotPhen.pApBdiff+ng1

##Plot frequency of recessive phenotype by generation number with different lines
#for initial allele frequencies (identical).
plotPhen.pApBsame<-ggplot(pApBvarySum_Same, aes(x = step, y = Phen,group = pAi)) + 
  ylab("Frequency of recessive
       phenotype")+xlab("Generation number")+
  geom_line(size = 0.5,aes(color = factor(pAi)))+
  coord_cartesian(ylim = c(0,1.0))+
  scale_y_continuous(breaks=seq(from = 0, to = 1.0, by = 0.1))+
  labs(colour = "Allele 
frequency")
plotPhen.pApBsame+ng1

#Dataset containing the first row from each simulation where the recessive phenotype is fixed
pApB_Phen.fixU <- datpApBVary[datpApBVary$Phen == 1.0,]
pApB_Phen.fixU <- pApB_Phen.fixU[,unique(pApB_Phen.fixU, by = c('pA.pB','sim'))]

#Dataset containing the first row from each simulation where the dominant phenotype is fixed
pApB_Phen.fixD <- datpApBVary[datpApBVary$Phen == 0.0,]
pApB_Phen.fixD <- pApB_Phen.fixD[,unique(pApB_Phen.fixD, by = c('pA.pB','sim'))]

#Datasets for number of simulations resulting in fixation of either phenotpype
pApB_fixU <- ddply(pApB_Phen.fixU,.(pA.pB),summarise,prop=(length(pA.pB)/1000))
pApB_fixD <- ddply(pApB_Phen.fixD,.(pA.pB),summarise,prop=(length(pA.pB)/1000))
pApB_fixD <- rbind(pApB_fixD,c('0.1;0.1',0.000),c('0.2;0.2',0.000),c('0.3;0.3',0.000),c('0.1;0.9',0.000),c('0.2;0.8',0.000))
pApB_fix <- merge(pApB_fixU,pApB_fixD,by.x ="pA.pB",by.y ="pA.pB")
setnames(pApB_fix,old=c("prop.x","prop.y"),new=c("FixU","FixD"))
pApB_fix$FixTot <- as.numeric(pApB_fix$FixU) + as.numeric(pApB_fix$FixD)

#Split pA.pB column into separate columns of allele frequencies and calculate difference in allele frequencies
pApB_fix <- tidyr::separate(pApB_fix, pA.pB, into = c("pA","pB"), ,remove = FALSE, sep = ";")
pApB_fix$Diff <- as.numeric(pApB_fix$pB) - as.numeric(pApB_fix$pA)

#Create two datasets: one for identical starting frequencies and one for different frequencies
pApB_fix_Diff <- pApB_fix[pApB_fix$Diff > 0,]
pApB_fix_Same <- pApB_fix[pApB_fix$Diff == 0,]

##Dataset showing mean time to fixation of recessive phenotype for each allele frequency combination
#Also mean frequency of both alleles at phenotype fixation
pApB_fixGen <- ddply(pApB_Phen.fixU,.(pA.pB),summarise,meanGen = mean(step),
                      mean.pA = mean(pA),mean.pB = mean(pB))

##Split pA.pB column from 'pApB_fixGen' dataset into separate columns 
#of allele frequencies and calculate difference in allele frequencies
pApB_fixGen <- tidyr::separate(pApB_fixGen, pA.pB, into = c("pA","pB"), ,remove = FALSE, sep = ";")
pApB_fixGen$Diff <- as.numeric(pApB_fixGen$pB) - as.numeric(pApB_fixGen$pA)

#Create two datasets: one for identical starting frequencies and one for different frequencies
pApB_fixGen_Diff <- pApB_fixGen[pApB_fixGen$Diff > 0,]
pApB_fixGen_Same <- pApB_fixGen[pApB_fixGen$Diff == 0,]

##Plot frequency of recessive phenotype by generation number with different lines
#for differences in initial allele frequencies.
plotFix.pApBdiff<-ggplot(pApB_fixGen_Diff, aes(x = factor(Diff), y = meanGen)) + 
  ylab("Mean generations
to fixation")+xlab("Difference in allele frequencies")+
  geom_bar(stat = "identity",color = "black", fill = "grey")+
  coord_cartesian(ylim = c(0,35))+
  scale_y_continuous(breaks=seq(from = 0, to = 35, by = 5))
plotFix.pApBdiff+ng1

##Plot frequency of recessive phenotype by generation number with different lines
#for differences in initial allele frequencies.
plotFix.pApBsame<-ggplot(pApB_fixGen_Same, aes(x = factor(pA), y = meanGen)) + 
  ylab("Mean generations
to fixation")+xlab("Starting allele frequencies")+
  geom_bar(stat = "identity",color = "black", fill = "grey")+
  coord_cartesian(ylim = c(0,45))+
  scale_y_continuous(breaks=seq(from = 0, to = 45, by = 5))
plotFix.pApBsame+ng1


## MODELS TO EXTRACT SLOPES (DIFFERENCE) ##

#Model for difference in allele frequencies = 0.2
pApBvarySum_Diff_0.2 <- subset(pApBvarySum_Diff, Diff == '0.2')
model_0.2 <- lm(Phen~step,data = pApBvarySum_Diff_0.2)
summary(model_0.2)
plot(Phen~step,data = pApBvarySum_Diff_0.2)
coefficients(model_0.2)[2]

#Model for difference in allele frequencies = 0.4
pApBvarySum_Diff_0.4 <- subset(pApBvarySum_Diff, Diff == '0.4')
model_0.4 <- lm(Phen~step,data = pApBvarySum_Diff_0.4)
summary(model_0.4)
plot(Phen~step,data = pApBvarySum_Diff_0.4)
coefficients(model_0.4)[2]

#Model for difference in allele frequencies = 0.6
pApBvarySum_Diff_0.6 <- subset(pApBvarySum_Diff, Diff == '0.6')
model_0.6 <- lm(Phen~step,data = pApBvarySum_Diff_0.6)
summary(model_0.6)
plot(Phen~step,data = pApBvarySum_Diff_0.6) # Looks non-linear
coefficients(model_0.6)[2]

#Model for difference in allele frequencies = 0.8
pApBvarySum_Diff_0.8 <- subset(pApBvarySum_Diff, Diff == '0.8')
model_0.8 <- lm(Phen~step,data = pApBvarySum_Diff_0.8)
summary(model_0.8)
plot(pB~step,data = pApBvarySum_Diff_0.8) # Looks non-linear
coefficients(model_0.8)[2]

