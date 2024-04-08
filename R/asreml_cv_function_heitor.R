################################################
#Cross-validation function for asreml GS models#
################################################
asreml_cv <- \(k,reps, data_file, G_matrix, response){
  acc_intervals <- matrix(NA,reps)
  h2_intervals <- matrix(NA,reps)
  for(s in 1:reps){
  n <- nrow(data_file)
  sel <- rep(1:k, length.out = n) 
  group <- sample(sel, n)
  ypred_cv <- matrix(data=NA, nrow=n, ncol=1)
  h2_cv <- matrix(data=NA, nrow=k, ncol=1)
  y <- response 
 
  for (g in 1:k)#the group that's going to be the test
  { 
   
    data_file$test <- y
    
    for (j in 1:n) {
      if(group[j] == g) { data_file$test[j]<-NA } #If an observation has the current "test group" card it is set to NA
    }

    modelGBLUPcv <- asreml(fixed = test ~ 1,
                           random = ~vm(genotype, G_matrix),
                           workspace = 128e06,
                           na.action = na.method(y = "include"),
                           data = data_file) 

    predGBLUPcv <- predict(modelGBLUPcv, classify="genotype", sed=T)$pvals
    predGBLUPcv <- predGBLUPcv[,2]
    
    (h2_GBLUPcv <- vpredict(modelGBLUPcv, h2~V1/(V1+V2)))
    (h2_cv[g] <- as.numeric(h2_GBLUPcv[1]))
    
    for (j in 1:n) {
      if(group[j] == g) { ypred_cv[j] <- predGBLUPcv[j] } #combines the predicted values of the test groups
    }
  }
  (PA <- cor(y, ypred_cv, method='pearson', use="complete.obs"))  
  
  acc_intervals[s] <- PA
  h2_intervals[s] <- mean(h2_cv)
  }
  return(list(accuracies = acc_intervals, h2 = h2_intervals))

}
