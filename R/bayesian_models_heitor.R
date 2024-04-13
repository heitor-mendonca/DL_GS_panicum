###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
#                                                                             #
#  Personal Project : use Deep Learning to perform GS in this article         #
###############################################################################
library(tidyverse)
library(BGLR)
#library(doParallel)
#library(magrittr)
#detectCores()
#-------------------------------------------------------------------------------
#Load data

load('./Data/molecular_data.RData')

geno <- X4[-26,] #removes genotype 29

fam_data <- readr::read_delim('./Data/blup_means.txt') %>%
  mutate(genotype = Offspring)%>% select(Mother,genotype) %>%
  mutate_all(as.factor) 

phenotype <- readr::read_rds('./Data/adj_mean_OM_heitor') %>% 
  inner_join(fam_data, by = c('genotype' ='genotype'))%>%
  mutate(family = Mother, resp = predicted.value) %>% 
  select(!c(Mother, predicted.value)) %>%  relocate(family) %>% droplevels()
#-------------------------------------------------------------------------------
# set some configs

nIter <- 6000; burnIn <- 700    # For Testing, actual training takes more

#-------------------------------------------------------------------------------
#Bayesian Ridge Regression

ETA <- list(list(X=geno, model="BRR"))
y <- phenotype$resp
gsBRR <- BGLR(y=y, ETA=ETA, nIter=nIter, burnIn=burnIn, saveAt = './Data/bayesian_files/BRR/')

#######################################################
# Calculating goodness-of-fit statistics (from salvador's script)

# Heritability (GS)
(Vary <- var(y, na.rm=TRUE))
(VarE <- gsBRR$varE)
(h2_GS <- 1 - VarE/Vary)

# Predictive Ability  corr(yadj,ghat)
predGS <- gsBRR$yHat
(PA <- cor(y, predGS, method='pearson', use="complete.obs")) # do NOT trust this, need to do CV

# Predictive Accuracy corr(greal,ghat)
# Need a vector with the greal
(ACC <- PA/sqrt(h2_GS))  # Approximation, obviously wrong










