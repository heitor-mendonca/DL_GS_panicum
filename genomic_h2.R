

###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
###############################################################################

rm(list=ls(all=T)) 
gc()

library(dplyr)
library(asreml)
library(asremlPlus)
library(AGHmatrix)


# Phenotipic data
means <- read.table(file="blup_means.txt", sep = "\t", dec = ".", header = TRUE)
head(means)


# Pedigree data
ped.data <- read.table("genealogy.txt", header=T, sep = "\t")
head(ped.data)
Amat <- Amatrix(data=ped.data, ploidy=4)
Amat[c(1:3,33:36), c(1:3,33:36)]
dim(Amat)

Amat_ind <- Amat[35:564,35:564]
dim(Amat_ind)                               # 530 x 530


# Genotipic data
Z <- diag(nrow(means))                      # Identity
dim(Z)                          

load("molecular_data.RData")
dim(X4)                                     # 530 offsprings and 41424 markers
Gmat <- Gmatrix(X4, method="VanRaden", ploidy=4, missingValue =NA)  
Gmat[1:5, 1:5]                   
dim(Gmat)                                   # 530 x 530


# Inverse
G <- 0.99*Gmat + 0.01*Amat_ind
Ginv <- solve(G)


## OM ####
data <- data.frame("ind" <- as.factor(means$Offspring), "OM" <- means$OM)
colnames(data) <- c("ind", "OM")

GBLUP_OM <- asreml(fixed = OM ~ 1, random = ~ vm(ind,G),
                   na.action = na.method(y = "include"), data = data)
GBLUP_OM <- update.asreml(GBLUP_OM)
infoCriteria.asreml(GBLUP_OM) 
summary(GBLUP_OM, all=T)$varcomp
summary(GBLUP_OM, all=T)$varcomp[1,1]/(summary(GBLUP_OM, all=T)$varcomp[1,1] + summary(GBLUP_OM, all=T)$varcomp[2,1])
# 0.5622

## IVD ####
data <- data.frame("ind" <- as.factor(means$Offspring), "IVD" <- means$IVD)
colnames(data) <- c("ind", "IVD")

GBLUP_IVD <- asreml(fixed = IVD ~ 1, random = ~ giv(ind),
                    ginverse = list(ind = Ginv),
                    na.method.Y = "include", data = data)
GBLUP_IVD <- update.asreml(GBLUP_IVD)
info.crit.asreml(GBLUP_IVD) 
summary(GBLUP_IVD, all=T)$varcomp
summary(GBLUP_IVD, all=T)$varcomp[1,2]/(summary(GBLUP_IVD, all=T)$varcomp[1,2] + summary(GBLUP_IVD, all=T)$varcomp[2,2])
# 0.2594

## CP ####
data <- data.frame("ind" <- as.factor(means$Offspring), "CP" <- means$CP)
colnames(data) <- c("ind", "CP")

GBLUP_CP <- asreml(fixed = CP ~ 1, random = ~ giv(ind),
                   ginverse = list(ind = Ginv),
                   na.method.Y = "include", data = data)
GBLUP_CP <- update.asreml(GBLUP_CP)
info.crit.asreml(GBLUP_CP) 
summary(GBLUP_CP, all=T)$varcomp
summary(GBLUP_CP, all=T)$varcomp[1,2]/(summary(GBLUP_CP, all=T)$varcomp[1,2] + summary(GBLUP_CP, all=T)$varcomp[2,2])
# 0.2656

## LDM ####
data <- data.frame("ind" <- as.factor(means$Offspring), "LDM" <- means$LDM)
colnames(data) <- c("ind", "LDM")

GBLUP_LDM <- asreml(fixed = LDM ~ 1, random = ~ giv(ind),
                    ginverse = list(ind = Ginv),
                    na.method.Y = "include", data = data)
GBLUP_LDM <- update.asreml(GBLUP_LDM)
info.crit.asreml(GBLUP_LDM) 
summary(GBLUP_LDM, all=T)$varcomp
summary(GBLUP_LDM, all=T)$varcomp[1,2]/(summary(GBLUP_LDM, all=T)$varcomp[1,2] + summary(GBLUP_LDM, all=T)$varcomp[2,2])
# 0.0423

## RC ####
data <- data.frame("ind" <- as.factor(means$Offspring), "RC" <- means$RC)
colnames(data) <- c("ind", "RC")

GBLUP_RC <- asreml(fixed = RC ~ 1, random = ~ giv(ind),
                   ginverse = list(ind = Ginv),
                   na.method.Y = "include", data = data)
GBLUP_RC <- update.asreml(GBLUP_RC)
info.crit.asreml(GBLUP_RC) 
summary(GBLUP_RC, all=T)$varcomp
summary(GBLUP_RC, all=T)$varcomp[1,2]/(summary(GBLUP_RC, all=T)$varcomp[1,2] + summary(GBLUP_RC, all=T)$varcomp[2,2])
# 0.3550

## PLB ####
data <- data.frame("ind" <- as.factor(means$Offspring), "PLB" <- means$PLB)
colnames(data) <- c("ind", "PLB")

GBLUP_PLB <- asreml(fixed = PLB ~ 1, random = ~ giv(ind),
                    ginverse = list(ind = Ginv),
                    na.method.Y = "include", data = data)
GBLUP_PLB <- update.asreml(GBLUP_PLB)
info.crit.asreml(GBLUP_PLB) 
summary(GBLUP_PLB, all=T)$varcomp
summary(GBLUP_PLB, all=T)$varcomp[1,2]/(summary(GBLUP_PLB, all=T)$varcomp[1,2] + summary(GBLUP_PLB, all=T)$varcomp[2,2])
# 0.3123


