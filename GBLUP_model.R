
###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
###############################################################################

rm(list=ls(all=T)) 
gc()

library(dplyr)
library(coda)
library(cvTools)                          
library(asreml)
library(asremlPlus)
library(nadiv)
library(AGHmatrix)


# Phenotipic data
means <- read.table(file="./blup_means.txt", sep = ",", dec = ".", header = TRUE)
head(means)
y.trait <- means$MOF


# Pedigree data
ped.data <- read.table("genealogy.txt", header=T, sep = "\t")
head(ped.data)
Amat <- Amatrix(data=ped.data, ploidy=4)
Amat[c(1:3,33:36), c(1:3,33:36)]
dim(Amat)

Amat_ind <- Amat[35:564,35:564]
dim(Amat_ind)                               # 530 x 530


# Genotipic data
Z <- diag(nrow(means))                      
dim(Z)                          

load("./molecular_data.RData")
dim(X4)                                     # 530 offsprings and 41424 markers
Gmat <- Gmatrix(X4, method="VanRaden", ploidy=4, missingValue=NA)  
Gmat[1:5, 1:5]                  
dim(Gmat)                                   # 530 x 530


# Inverse
G <- 0.99*Gmat + 0.01*Amat_ind
Ginv <- solve(G)


##Cross-validation sets:####
n_id = length(means$Offspring)               
n_model = 1
s = 1000

y_pred_list = NULL
y_pred_list = c(y_pred_list, 1)
y_pred_G = matrix(NA,n_id,s)


for(k in 1:s){
  ## Sets: ####  
  seed_value = as.integer(paste("2001",k,sep=""))     # Fix the seed to each iterations
  set.seed(seed_value)                      # Setting the seed
  cv = cvFolds(n_id, K = 5)                 # Obtaining cross-validation subsets
  cv_vals = as.matrix(unlist(cv["subsets"][1]))       # Indexing subsets
  
  ID_names = means$Offspring                # Getting the individuals names
  
  out_set1 = cv_vals[unlist(cv['which'])==1,1]        # Indexing subset1
  out_set2 = cv_vals[unlist(cv['which'])==2,1]        # Indexing subset2
  out_set3 = cv_vals[unlist(cv['which'])==3,1]        # Indexing subset3  
  out_set4 = cv_vals[unlist(cv['which'])==4,1]        # Indexing subset4
  out_set5 = cv_vals[unlist(cv['which'])==5,1]        # Indexing subset5
  
  out_sets = list(out_set1, out_set2, out_set3, out_set4, out_set5)
  
  set1 = c(out_set2,out_set3,out_set4,out_set5)       # out_set1 out
  set1 = set1[order(set1,decreasing=F)]
  set2 = c(out_set1,out_set3,out_set4,out_set5)       # out_set2 out
  set2 = set2[order(set2,decreasing=F)]
  set3 = c(out_set1,out_set2,out_set4,out_set5)       # out_set3 out
  set3 = set3[order(set3,decreasing=F)]
  set4 = c(out_set1,out_set2,out_set3,out_set5)       # out_set4 out
  set4 = set4[order(set4,decreasing=F)]
  set5 = c(out_set1,out_set2,out_set3,out_set4)       # out_set5 out
  set5 = set5[order(set5,decreasing=F)]
  
  # Evaluating with there is some ID common across out_set
  all(!c(out_set1%in%out_set2,out_set1%in%out_set3,out_set1%in%out_set4,out_set1%in%out_set5,out_set2%in%out_set3,out_set2%in%out_set4,out_set2%in%out_set5,out_set3%in%out_set4,out_set3%in%out_set5,out_set4%in%out_set5)) #OK!
  
  sets = list(set1=set1,set2=set2,set3=set3,set4=set4,set5=set5)
  out_sets = list(out_set1=out_set1,out_set2=out_set2,out_set3=out_set3,out_set4=out_set4,out_set5=out_set5)
  rm(out_set1,out_set2,out_set3,out_set4,out_set5,set1,set2,set3,set4,set5)
  
  ## G-BLUP: ####
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set1]=y.trait[sets$set1]
  data <- data.frame("ind" <- as.factor(means$Offspring), "trait" <- y_obs)
  colnames(data) <- c("ind", "trait")
  GBLUP_set1 <- asreml(fixed = trait ~ 1, random = ~ giv(ind),
                       ginverse = list(ind = Ginv),
                       na.method.Y = "include", data = data)
  GBLUP_set1 <- update.asreml(GBLUP_set1)
  pred_set1 <- predict(GBLUP_set1, classify="ind", sed=T)$predictions$pvals
  y_pred_G[out_sets$out_set1, k] <- pred_set1[out_sets$out_set1, 2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set2]=y.trait[sets$set2]
  data <- data.frame("ind" <- as.factor(means$Offspring), "trait" <- y_obs)
  colnames(data) <- c("ind", "trait")
  GBLUP_set2 <- asreml(fixed = trait ~ 1, random = ~ giv(ind),
                         ginverse = list(ind = Ginv),
                         na.method.Y = "include", data = data)
  GBLUP_set2 <- update.asreml(GBLUP_set2)
  pred_set2 <- predict(GBLUP_set2, classify="ind", sed=T)$predictions$pvals
  y_pred_G[out_sets$out_set2, k] <- pred_set2[out_sets$out_set2, 2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set3]=y.trait[sets$set3]
  data <- data.frame("ind" <- as.factor(means$Offspring), "trait" <- y_obs)
  colnames(data) <- c("ind", "trait")
  GBLUP_set3 <- asreml(fixed = trait ~ 1, random = ~ giv(ind),
                       ginverse = list(ind = Ginv),
                       na.method.Y = "include", data = data)
  GBLUP_set3 <- update.asreml(GBLUP_set3)
  pred_set3 <- predict(GBLUP_set3, classify="ind", sed=T)$predictions$pvals
  y_pred_G[out_sets$out_set3, k] <- pred_set3[out_sets$out_set3, 2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set4]=y.trait[sets$set4]
  data <- data.frame("ind" <- as.factor(means$Offspring), "trait" <- y_obs)
  colnames(data) <- c("ind", "trait")
  GBLUP_set4 <- asreml(fixed = trait ~ 1, random = ~ giv(ind),
                         ginverse = list(ind = Ginv),
                         na.method.Y = "include", data = data)
  GBLUP_set4 <- update.asreml(GBLUP_set4)
  pred_set4 <- predict(GBLUP_set4, classify="ind", sed=T)$predictions$pvals
  y_pred_G[out_sets$out_set4, k] <- pred_set4[out_sets$out_set4, 2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set5]=y.trait[sets$set5]
  data <- data.frame("ind" <- as.factor(means$Offspring), "trait" <- y_obs)
  colnames(data) <- c("ind", "trait")
  GBLUP_set5 <- asreml(fixed = trait ~ 1, random = ~ giv(ind),
                         ginverse = list(ind = Ginv),
                         na.method.Y = "include", data = data)
  GBLUP_set5 <- update.asreml(GBLUP_set5)
  pred_set5 <- predict(GBLUP_set5, classify="ind", sed=T)$predictions$pvals
  y_pred_G[out_sets$out_set5, k] <- pred_set5[out_sets$out_set5, 2]
  
  y_pred_list = c(y_pred_list,c(y_pred_G[,k]))
}

all(is.na(y_pred_list))


for (i in 1:s) {
  if(i==1) {
    accuracy_intervals = matrix(NA,s,n_model)
    press = matrix(NA,s,n_model)
  }
  accuracy_intervals[i,1] = cor(y_pred_G[,i], y.trait)
  press[i,1] = sum((y.trait - y_pred_G[,i])^2)
}

dim(accuracy_intervals)
dim(press)


colnames(accuracy_intervals) <- c("GBLUP")
accuracy_intervals_mcmc = as.mcmc(accuracy_intervals)
HPDinterval(accuracy_intervals_mcmc, prob=0.95)
colMeans(accuracy_intervals) 

colnames(press) <- c("GBLUP")
colMeans(press)

