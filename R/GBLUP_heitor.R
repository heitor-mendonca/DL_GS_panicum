###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
#                                                                             #
#  Personal Project : use Deep Learning to perform GS in this article         #
###############################################################################
library(tidyverse)
library(asreml)
library(asremlPlus)
library(ASRgenomics)
library(coda)
source("./R/asreml_cv_function_heitor.R")
#Loading required data


G_blend <- readr::read_rds('./Data/G_blend_heitor') #99% G mat
G_clean <- readr::read_rds('./Data/G_clean_heitor')
G_align <- readr::read_rds('./Data/G_aligned_heitor')
H <-       readr::read_rds('./Data/H_heitor')
fam_data <- readr::read_delim('./Data/blup_means.txt') %>%
  mutate(genotype = Offspring)%>% select(Mother,genotype) %>%
  mutate_all(as.factor) 
phenotype <- readr::read_rds('./Data/adj_mean_OM_heitor') %>% 
  inner_join(fam_data, by = c('genotype' ='genotype'))%>%
  mutate(family = Mother, resp = predicted.value) %>% 
  select(!c(Mother, predicted.value)) %>%  relocate(family) %>% droplevels()


str(phenotype)
str(fam_data)

#-------------------------------------------------------------------------------
#
check <- match.kinship2pheno(K=G_mat, pheno.data=phenotype,
                             indiv='genotype', clean=FALSE, mism=TRUE)
#-------------------------------------------------------------------------------
#run CV function for GBLUP with blended matrix

cv_gblup <- asreml_cv(k = 5,reps= 5000, data_file = phenotype, G_matrix = G_mat,
                      response = phenotype$resp)

gblup_acc <- cv_gblup[[1]]
gblup_h2  <- cv_gblup[[2]]

accuracy_intervals_mcmc = as.mcmc(gblup_acc)
HPDinterval(accuracy_intervals_mcmc, prob=0.95)
colMeans(accuracy_intervals_mcmc) 

#lara's hpdi 38-44  0.4168932 
#mine 39-45  0.4305951 (with 100 iter) 0.3942776-0.4589232  0.4268566
#-------------------------------------------------------------------------------
# GBLUP with raw G
cv_gblup_clean <- asreml_cv(k = 5,reps= 100, data_file = phenotype,
                            G_matrix = G_clean,response = phenotype$resp)

gblup_clean_acc <- cv_gblup_clean[[1]]
gblup_clean_h2  <- cv_gblup_clean[[2]]


accuracy_intervals_mcmc = as.mcmc(gblup_clean_acc)
HPDinterval(accuracy_intervals_mcmc, prob=0.95)
colMeans(accuracy_intervals_mcmc) 
