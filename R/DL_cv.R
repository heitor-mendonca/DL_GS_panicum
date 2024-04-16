units_M <- 200 #neurons per hidden layer
epochs_M <- 30
reps <- 10

acc_intervals <- matrix(NA,reps)
for(s in 1:reps){
ypred_cv <- matrix(data=NA, nrow=n, ncol=1)
y <- phenotype$resp
n <- nrow(chol_blend)
sel <- rep(1:5, length.out = n) 
group <- sample(sel, n)


for(i in 1:5){
X_tr <- chol_blend[!group %in% i,]
X_ts <- chol_blend[group %in% i,]
y_tr <- scale(y[!group %in% i])

model <- keras_model_sequential() %>%
  layer_dense(
    units = units_M,
    activation = 'relu',
    input_shape = c(dim(X_tr)[2]))%>%
  layer_dropout(rate=.2)%>%#input
  layer_dense(units = units_M, activation = 'relu')%>%
  layer_dropout(rate=.2)%>%#hidden 1
  layer_dense(units = units_M, activation = 'relu')%>%
  layer_dropout(rate=.2)%>%#hidden 2
  layer_dense(units = 1)#output


model %>% compile(loss = 'mean_squared_error',
                  optimizer = optimizer_adam(),
                  metrics = c('mean_squared_error'))  

summary(model)

history <- model %>% fit(
  X_tr,y_tr, epochs = epochs_M, batch_size = 64,validation_split = .2
)

y_p <- model %>% predict(X_ts)
ypred_cv[group%in%i] <- y_p
}
(PA <- cor(y, ypred_cv, method='pearson', use="complete.obs"))  
acc_intervals[s] <- PA
}


accuracy_intervals_mcmc = as.mcmc(acc_intervals)
HPDinterval(accuracy_intervals_mcmc)
mean(acc_intervals)
cor(y, ypred_cv)

