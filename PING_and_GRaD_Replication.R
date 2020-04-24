#Explanation: We discovered a series of SNPs that were associated with decoding
#performance, an early predictor of reading ability. When we went to replicate
#our results, we were only able to replicate if we age matched the cohort. Here
#I am investigating if there is an interaction between top SNP and age


#Andrew K. Adans

#Load a bunch of needed packages
library(tidyverse)
library(data.table)
library(stargazer)
library(sjPlot)

#set the local directory for the project
baseDir <- "/Volumes/NHLP_WGS/NHLP_Version2_Analysis/Replications/"

colRename <- function(df, x, name) {
  #rename the first column in a dataframe and return that frame 
  names(df)[x] <- c(name)
  return(df)
}

#GRaD Replication
#PLINK format files provide a family id column, we don't need it, drop before combining
gradChr19Geno <- fread(file = paste(baseDir, 'GRaD_Replication/Chr19_Top_Peak.raw', sep = ''),
                       sep = ' ', header = TRUE)
gradCovs <- fread(file = paste(baseDir, "GRaD_Replication/FullGRaD_Covariate_Inattention_EV50.cov", sep = ''), 
                  na.strings = '-9', sep = '\t', header = TRUE)
fullGradPheno <- fread(file = paste(baseDir, 'GRaD_Replication/GradReadLatent.txt', sep = ''), 
                       na.strings = '-9', sep = ' ', header = TRUE)
ageMatchGradPheno <- fread(file = paste(baseDir, "GRad_Replication/GradReadLatentYoung.txt", sep = ''),
                          na.strings = '-9', sep = ' ', header = TRUE)
gradCovs$SES_CC <- recode(gradCovs$SES_CC, '1' = 0, '2' = 1)
gradChr19Geno$SEX <- recode(gradChr19Geno$SEX, '1' = 0, '2' = 1)

#Merge full GRaD for moderation models
fullGrad <- list(fullGradPheno, gradCovs, gradChr19Geno)
fullGrad <- lapply(fullGrad, colRename, x = 2, name = 'gradID')
fullGradJoin <- plyr::join_all(fullGrad, by = 'gradID', type = 'inner')
fullGradJoin$age <- as.numeric(fullGradJoin$age)

#Merge up the age matched version
ageMatchedGrad <- list(ageMatchGradPheno, gradCovs, gradChr19Geno)
ageMatchedGrad <- lapply(ageMatchedGrad, colRename, x= 2, name = 'gradID')
ageMatchedGradJoin <- plyr::join_all(ageMatchedGrad, by = 'gradID', type = 'inner')

#Replicate NHLP top SNP in GRaD
ageMatchedGradJoin$ZGradReadLatentYoungCC <- recode(ageMatchedGradJoin$ZGradReadLatentYoungCC, '1' = 0, '2' = 1)

gradRep <- glm(ZGradReadLatentYoungCC ~ rs2599553_A + TestAgeTotalMonths + SEX + SES_CC + 
                 EV1 + EV2 + EV3 + EV4 + EV5 + EV6 + EV7 + EV8 + EV9 + EV10, 
               data = ageMatchedGradJoin, family = binomial(link = "logit"))
summary(gradRep)

#Does the quantitative version give significant results for the full GRaD cohort
gradLin <- lm(GradReadLatent ~ rs2599553_A + age + SEX + SES_CC +
                EV1 + EV2 + EV3 + EV4 + EV5 + EV6 + EV7 + EV8 + EV9 + EV10,
               data = fullGradJoin)
summary(gradLin)

gradMod <- lm(GradReadLatent ~ rs2599553_A*age + SEX + SES_CC +
                EV1 + EV2 + EV3 + EV4 + EV5 + EV6 + EV7 + EV8 + EV9 + EV10, 
               data = fullGradJoin)
summary(gradMod)
stargazer(gradLin, gradMod, type = 'text',
          column.labels = c('No Interation', 'Age/SNP Interaction'),
          dep.var.labels = 'GRaD Latent Reading Variable',
          model.numbers = FALSE,
          ci = TRUE,
          intercept.bottom = FALSE,
          star.char = c('.', '*', '**'),
          single.row = TRUE,
          out = paste(baseDir, 'GRaD_Replication/GRaD_all_ages.txt', sep = ''))

#Plot the main effects by quartile
plot_model(gradMod, type = "int", terms = c("rs2599553_A [0,1,2]", "age") , mdrt.values = "quart")

