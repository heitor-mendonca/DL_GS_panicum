
###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
###############################################################################

rm(list=ls(all=T)) 
gc()

library(foreach)
library(doMC)
registerDoMC(32)               # be aware with this command line -> it is using 32 cores
library(BGLR)
library(dplyr)
library(coda)
library(cvTools)
library(AGHmatrix)


# Phenotipic data
means <- read.table(file="./blup_means.txt", sep = ",", dec = ".", header = TRUE)
head(means)
y.trait <- means$MOF


# Genotipic data
Z <- diag(nrow(means))                      
dim(Z)                          

load("./molecular_data.RData")
dim(X4)                                     # 530 offsprings and 41424 markers

nIter <- 20000 
burnIn <- 2000


##Cross-validation sets:####
ETA.BRR <- list(X=list(X=X4, model='BRR'))           # BRR -> Gaussian Prior
ETA.BA <- list(X=list(X=X4, model='BayesA'))         # Bayes A -> Scaled-t Prior 
ETA.BB <- list(X=list(X=X4, model='BayesB'))         # Bayes B -> Point of mass at zero + scaled-t slab
ETA.BC <- list(X=list(X=X4, model='BayesC'))         # Bayes C -> Point of mass at zero + Gaussian slab
ETA.BL <- list(X=list(X=X4, model='BL'))             # BLASSO -> Laplace prior 

n_id = length(means$Offspring)              
n_model = 5
s = 100

y_pred_list = NULL
y_pred_list = c(y_pred_list, 1)
y_pred.BRR = matrix(NA,n_id,s)
y_pred.BA = matrix(NA,n_id,s)
y_pred.BB = matrix(NA,n_id,s)
y_pred.BC = matrix(NA,n_id,s)
y_pred.BL = matrix(NA,n_id,s)


list_for = foreach(k=1:s, .export="y_pred_list") %dopar% {
  ## Sets: ####  
  seed_value = as.integer(paste("2001",k,sep=""))     # Fix the seed to each iterations
  set.seed(seed_value)                      
  cv = cvFolds(n_id, K = 5)                 # Obtaining cross-validation subsets
  cv_vals = as.matrix(unlist(cv["subsets"][1]))       # Indexing subsets
  
  ID_names = means$Offspring                
  
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
  
  ## BRR: ####
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set1]=y.trait[sets$set1]
  ans.brr1 <- BGLR(y=y_obs, ETA=ETA.BRR, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BRR[out_sets$out_set1,k] <- ans.brr1$yHat[out_sets$out_set1]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set2]=y.trait[sets$set2]
  ans.brr2 <- BGLR(y=y_obs, ETA=ETA.BRR, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BRR[out_sets$out_set2,k] <- ans.brr2$yHat[out_sets$out_set2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set3]=y.trait[sets$set3]
  ans.brr3 <- BGLR(y=y_obs, ETA=ETA.BRR, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BRR[out_sets$out_set3,k] <- ans.brr3$yHat[out_sets$out_set3]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set4]=y.trait[sets$set4]
  ans.brr4 <- BGLR(y=y_obs, ETA=ETA.BRR, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BRR[out_sets$out_set4,k] <- ans.brr4$yHat[out_sets$out_set4]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set5]=y.trait[sets$set5]
  ans.brr5 <- BGLR(y=y_obs, ETA=ETA.BRR, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BRR[out_sets$out_set5,k] <- ans.brr5$yHat[out_sets$out_set5]

  y_pred_list = c(y_pred_list,c(y_pred.BRR[,k]))
  
  ## Bayes A: ####
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set1]=y.trait[sets$set1]
  ans.ba1 <- BGLR(y=y_obs, ETA=ETA.BA, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BA[out_sets$out_set1,k] <- ans.ba1$yHat[out_sets$out_set1]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set2]=y.trait[sets$set2]
  ans.ba2 <- BGLR(y=y_obs, ETA=ETA.BA, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BA[out_sets$out_set2,k] <- ans.ba2$yHat[out_sets$out_set2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set3]=y.trait[sets$set3]
  ans.ba3 <- BGLR(y=y_obs, ETA=ETA.BA, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BA[out_sets$out_set3,k] <- ans.ba3$yHat[out_sets$out_set3]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set4]=y.trait[sets$set4]
  ans.ba4 <- BGLR(y=y_obs, ETA=ETA.BA, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BA[out_sets$out_set4,k] <- ans.ba4$yHat[out_sets$out_set4]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set5]=y.trait[sets$set5]
  ans.ba5 <- BGLR(y=y_obs, ETA=ETA.BA, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BA[out_sets$out_set5,k] <- ans.ba5$yHat[out_sets$out_set5]

  y_pred_list = c(y_pred_list,c(y_pred.BA[,k]))
  
  ## Bayes B: ####
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set1]=y.trait[sets$set1]
  ans.bb1 <- BGLR(y=y_obs, ETA=ETA.BB, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BB[out_sets$out_set1,k] <- ans.bb1$yHat[out_sets$out_set1]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set2]=y.trait[sets$set2]
  ans.bb2 <- BGLR(y=y_obs, ETA=ETA.BB, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BB[out_sets$out_set2,k] <- ans.bb2$yHat[out_sets$out_set2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set3]=y.trait[sets$set3]
  ans.bb3 <- BGLR(y=y_obs, ETA=ETA.BB, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BB[out_sets$out_set3,k] <- ans.bb3$yHat[out_sets$out_set3]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set4]=y.trait[sets$set4]
  ans.bb4 <- BGLR(y=y_obs, ETA=ETA.BB, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BB[out_sets$out_set4,k] <- ans.bb4$yHat[out_sets$out_set4]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set5]=y.trait[sets$set5]
  ans.bb5 <- BGLR(y=y_obs, ETA=ETA.BB, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BB[out_sets$out_set5,k] <- ans.bb5$yHat[out_sets$out_set5]

  y_pred_list = c(y_pred_list,c(y_pred.BB[,k]))
  
  ## Bayes C: ####
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set1]=y.trait[sets$set1]
  ans.bc1 <- BGLR(y=y_obs, ETA=ETA.BC, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BC[out_sets$out_set1,k] <- ans.bc1$yHat[out_sets$out_set1]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set2]=y.trait[sets$set2]
  ans.bc2 <- BGLR(y=y_obs, ETA=ETA.BC, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BC[out_sets$out_set2,k] <- ans.bc2$yHat[out_sets$out_set2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set3]=y.trait[sets$set3]
  ans.bc3 <- BGLR(y=y_obs, ETA=ETA.BC, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BC[out_sets$out_set3,k] <- ans.bc3$yHat[out_sets$out_set3]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set4]=y.trait[sets$set4]
  ans.bc4 <- BGLR(y=y_obs, ETA=ETA.BC, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BC[out_sets$out_set4,k] <- ans.bc4$yHat[out_sets$out_set4]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set5]=y.trait[sets$set5]
  ans.bc5 <- BGLR(y=y_obs, ETA=ETA.BC, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BC[out_sets$out_set5,k] <- ans.bc5$yHat[out_sets$out_set5]

  y_pred_list = c(y_pred_list,c(y_pred.BC[,k]))
  
  ## BLASSO: ####
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set1]=y.trait[sets$set1]
  ans.bl1 <- BGLR(y=y_obs, ETA=ETA.BL, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BL[out_sets$out_set1,k] <- ans.bl1$yHat[out_sets$out_set1]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set2]=y.trait[sets$set2]
  ans.bl2 <- BGLR(y=y_obs, ETA=ETA.BL, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BL[out_sets$out_set2,k] <- ans.bl2$yHat[out_sets$out_set2]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set3]=y.trait[sets$set3]
  ans.bl3 <- BGLR(y=y_obs, ETA=ETA.BL, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BL[out_sets$out_set3,k] <- ans.bl3$yHat[out_sets$out_set3]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set4]=y.trait[sets$set4]
  ans.bl4 <- BGLR(y=y_obs, ETA=ETA.BL, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BL[out_sets$out_set4,k] <- ans.bl4$yHat[out_sets$out_set4]
  
  y_obs=matrix(NA,n_id,1)
  y_obs[sets$set5]=y.trait[sets$set5]
  ans.bl5 <- BGLR(y=y_obs, ETA=ETA.BL, nIter=nIter, burnIn=burnIn, thin=5, verbose = F)
  y_pred.BL[out_sets$out_set5,k] <- ans.bl5$yHat[out_sets$out_set5]
  
  y_pred_list = c(y_pred_list,c(y_pred.BL[,k]))
  
}

all(is.na(y_pred_list))

y_pred_matrix = matrix(NA,(n_id*n_model),s)
names = paste("G", rep(1:n_id,times=n_model),"_M",rep(1:n_model,each=n_id),sep="")
rownames(y_pred_matrix) <- names
colnames(y_pred_matrix) <- paste("s",1:s,sep="")

count = 1
for (i in 1:s) {
  y_pred_matrix[,i] = list_for[[i]][-1]
}

y_pred.BRR_final = y_pred_matrix[grepl(c("M1"),rownames(y_pred_matrix)),]
y_pred.BA_final = y_pred_matrix[grepl(c("M2"),rownames(y_pred_matrix)),]
y_pred.BB_final = y_pred_matrix[grepl(c("M3"),rownames(y_pred_matrix)),]
y_pred.BC_final = y_pred_matrix[grepl(c("M4"),rownames(y_pred_matrix)),]
y_pred.BL_final = y_pred_matrix[grepl(c("M5"),rownames(y_pred_matrix)),]


for (i in 1:s) {
  if(i==1) {
    accuracy_intervals = matrix(NA,s,n_model)
    press = matrix(NA,s,n_model)
  }
  accuracy_intervals[i,1] = cor(y_pred.BRR_final[,i], y.trait)
  accuracy_intervals[i,2] = cor(y_pred.BA_final[,i], y.trait)
  accuracy_intervals[i,3] = cor(y_pred.BB_final[,i], y.trait)
  accuracy_intervals[i,4] = cor(y_pred.BC_final[,i], y.trait)
  accuracy_intervals[i,5] = cor(y_pred.BL_final[,i], y.trait)

  press[i,1] = sum((y.trait - y_pred.BRR_final[,i])^2)
  press[i,2] = sum((y.trait - y_pred.BA_final[,i])^2)
  press[i,3] = sum((y.trait - y_pred.BB_final[,i])^2)
  press[i,4] = sum((y.trait - y_pred.BC_final[,i])^2)
  press[i,5] = sum((y.trait - y_pred.BL_final[,i])^2)
}
dim(accuracy_intervals)

colnames(accuracy_intervals) <- c("BRR","BayesA", "BayesB", "BayesC", "BLASSO")
accuracy_intervals_mcmc = as.mcmc(accuracy_intervals)
HPDinterval(accuracy_intervals_mcmc, prob=0.95)
colMeans(accuracy_intervals)

colnames(press) <- c("BRR","BayesA", "BayesB", "BayesC", "BLASSO")
colMeans(press)


