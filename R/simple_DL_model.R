library(tidyverse)
library(keras)
#------------------------------------------------------------------------
#import data

G_blend <- readr::read_rds('./Data/G_blend_heitor') #99% G mat
#G_clean <- readr::read_rds('./Data/G_clean_heitor')
#G_align <- readr::read_rds('./Data/G_aligned_heitor')
#H <-       readr::read_rds('./Data/H_heitor')
fam_data <- readr::read_delim('./Data/blup_means.txt') %>%
  mutate(genotype = Offspring)%>% select(Mother,genotype) %>%
  mutate_all(as.factor) 
phenotype <- readr::read_rds('./Data/adj_mean_OM_heitor') %>% 
  inner_join(fam_data, by = c('genotype' ='genotype'))%>%
  mutate(family = Mother, resp = predicted.value) %>% 
  select(!c(Mother, predicted.value)) %>%  relocate(family) %>% droplevels()#%>%
  #arrange(as.integer(genotype)) 
#-------------------------------------------------------------------------------
#Prepare data, create training and test sets
chol_blend <- base::chol(G_blend) %>% t()

str(phenotype)
#Z1G <- model.matrix(~-1+as.factor(phenotype$genotype))
#Z1Gt <- Z1G %*% chol_blend
y <- phenotype$resp
#X <- Z1G
n <- dim(G_blend)[1]
Post_trn <- sample(1:n,round(n*0.8))
X_tr <- chol_blend[Post_trn,]
X_ts <- chol_blend[-Post_trn,]
y_tr <- scale(y[Post_trn])
y_ts<- scale(y[-Post_trn])
Mean_trn <- mean(y[Post_trn])
SD_trn <- sd((y[Post_trn]))
#y_ts<- (y[-Post_trn] - Mean_trn)/SD_trn 
#-------------------------------------------------------------------------------
#specify net architecture, simple MLP

units_M <- 220 #neurons per hidden layer
epochs_M <- 20


model <- keras_model_sequential() %>%
  layer_dense(
    units = units_M,
    activation = 'relu',
    input_shape = c(dim(X_tr)[2]))%>%
  layer_dropout(rate=.3)%>%#input
  layer_dense(units = units_M, activation = 'relu')%>%
  layer_dropout(rate=.3)%>%#hidden 1
  layer_dense(units = units_M, activation = 'relu')%>%
  layer_dropout(rate=.3)%>%#hidden 2
  layer_dense(units = 1)#output


model %>% compile(loss = 'mean_squared_error',
                  optimizer = optimizer_adam(),
                  metrics = c('mean_squared_error'))  

summary(model)

history <- model %>% fit(
  X_tr,y_tr, epochs = epochs_M, batch_size = 50,validation_split = .2
)

#----------------------------------------------------------------
#Assess model performance
pf <- model %>% evaluate(x=X_ts, y = y_ts)

y_p <- model %>% predict(X_ts)
y_p = y_p * SD_trn + Mean_trn
y_ts = y_ts
y_ts = y_ts *SD_trn + Mean_trn


Y_all_tst <- data.frame(cbind(y_ts,y_p))
cor(y_ts,y_p)
plot(y_ts,y_p)
#---------------------------------------------------------------
# Trying with raw marker data, I'll use Lara's adjusted means
rm(list=ls())

load('./Data/molecular_data.RData')
geno <- X4

phenotype <- readr::read_delim('./Data/blup_means.txt') %>%
  mutate(genotype = Offspring)%>% select(Mother,genotype, OM) 

#-------------------------------------------------------------------------------
#Prepare data, create training and test sets

y <- phenotype$OM
n <- dim(geno)[1]
Post_trn <- sample(1:n,round(n*0.8))
X_tr <- geno[Post_trn,]
X_ts <- geno[-Post_trn,]
y_tr <- scale(y[Post_trn])
y_ts<- scale(y[-Post_trn])
Mean_trn <- mean(y[Post_trn])
SD_trn <- sd((y[Post_trn]))
#y_ts<- (y[-Post_trn] - Mean_trn)/SD_trn 
#-------------------------------------------------------------------------------
#specify net architecture, simple MLP

units_M <- 220 #neurons per hidden layer
epochs_M <- 20


model <- keras_model_sequential() %>%
  layer_dense(
    units = units_M,
    activation = 'relu',
    input_shape = c(dim(X_tr)[2]))%>%
  layer_dropout(rate=.3)%>%#input
  layer_dense(units = units_M, activation = 'relu')%>%
  layer_dropout(rate=.3)%>%#hidden 1
  layer_dense(units = units_M, activation = 'relu')%>%
  layer_dropout(rate=.3)%>%#hidden 2
  layer_dense(units = 1)#output


model %>% compile(loss = 'mean_squared_error',
                  optimizer = optimizer_adam(),
                  metrics = c('mean_squared_error'))  

summary(model)

history <- model %>% fit(
  X_tr,y_tr, epochs = epochs_M, batch_size = 50,validation_split = .2
)

#----------------------------------------------------------------
#Assess model performance
pf <- model %>% evaluate(x=X_ts, y = y_ts)

y_p <- model %>% predict(X_ts)
y_p = y_p * SD_trn + Mean_trn
y_ts = y_ts
y_ts = y_ts *SD_trn + Mean_trn


Y_all_tst <- data.frame(cbind(y_ts,y_p))
cor(y_ts,y_p)
plot(y_ts,y_p)
