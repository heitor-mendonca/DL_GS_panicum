###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
#                                                                             #
#  Personal Project : use Deep Learning to perform GS in this article         #
###############################################################################
library(AGHmatrix)
library(tidyverse)
library(ASRgenomics)
library(asreml)

#load molecular_data
load('./Data/molecular_data.RData')
geno <- X4


#--------------------------------------------------------------------------
#Obtaining G 

G <- AGHmatrix::Gmatrix(geno, ploidy = 4,missingValue =NA)

check_G <- kinship.diagnostics(G)
check_G$list.diagonal
check_G$plot.diag
#---------------------------------------------------------------------------
#Obtaining A

pedigree <- read.table('./Data/genealogy.txt', header = T, sep = '\t') #don't use the readr package here
A <- AGHmatrix::Amatrix(pedigree, ploidy = 4)
#----------------------------------------------------------------------------
#comparing G and A

G2A <- match.G2A(A=A, G=G, clean=TRUE, ord=TRUE, mism=TRUE, RMdiff=TRUE)

G2A$plotG2A

#misaligned

#-----------------------------------------------------------------------------
#manipulating G


#Aligned G

Gclean<- G2A$Gclean #intersect of G and A
Aclean<- G2A$Aclean

G_aligned <- G.tuneup(Gclean,Aclean, align = T) %>% `[[`('Gb')
G_aligned_inv <- G.inverse(G_aligned, sparseform = T)$Ginv
kinship.diagnostics(G_aligned) #less extreme diagonal values


#Blended G (equivalent to what Lara did in her script)

G_blend <- G.tuneup(Gclean, Aclean, blend = T, pblend = .01)%>% `[[`('Gb')
dim(G_blend)

G_blend_inv <- G.inverse(G_blend, sparseform = T)$Ginv
kinship.diagnostics(G_blend) #more extreme diagonal values than aligned, but mean around 1

#------------------------------------------------------------------------------
#Obtaining H matrix
Ginv <- G.inverse(G = G_blend, sparseform = FALSE)$Ginv
Ginv[1:5, 1:5] 
Hinv.sparse <- H.inverse(A = A, G = Ginv, lambda = 0.90, sparseform = TRUE)
head(Hinv.sparse)

H <- H.matrix(A = A, Ginv = Ginv, lambda = 0.90, sparseform = FALSE)
dim(H)



A[25:30, 25:30]
H[25:30, 25:30]
A[44:49, 44:49]
H[44:49, 44:49]
A[505:510, 505:510]
H[505:510, 505:510]
#------------------------------------------------------------------------------
#obtaining genomic heritability 

#import means, intersect with G_align

means <- readr::read_delim('./Data/adjusted_means_heitor.txt') %>% 
  mutate(genotype = as.character(genotype)) %>%
  filter(genotype %in% colnames(G_aligned)) %>% 
  mutate(genotype = as.factor(genotype))

str(means)

#fit model with aligned matrix

mod_h2 <- asreml(fixed = predicted.value ~ 1,
                 random = ~ vm(genotype, G_aligned_inv),
                 data = means)
summary(mod_h2)$varcomp
vpredict(mod_h2, H2 ~ V1/(V1+V2)) #narrow sense h2 0.5962427; broad sense h2 0.6930863

#fit model with blended matrix

mod_h2 <- asreml(fixed = predicted.value ~ 1,
                 random = ~ vm(genotype, G_blend_inv),
                 data = means)


summary(mod_h2)$varcomp

vpredict(mod_h2, H2 ~ V1/(V1+V2)) #narrow sense h2 0.5812985 

#fit model with H matrix

#needs more complete phenotype
means_H <- readr::read_delim('./Data/adjusted_means_heitor.txt') %>% 
  mutate(genotype = as.character(genotype))%>%
  filter(genotype %in% colnames(H)) %>% 
  mutate(genotype = as.factor(genotype))


mod_h2 <- asreml(fixed = predicted.value ~ 1,
                 random = ~ vm(genotype, Hinv.sparse),
                 data = means_H)


summary(mod_h2)$varcomp

vpredict(mod_h2, H2 ~ V1/(V1+V2)) #narrow sense h2 0.6390203 
#------------------------------------------------------------------------------
#saving matrices and phentoypes for GS

readr::write_rds(G_aligned, './Data/G_aligned_heitor')
readr::write_rds(Gclean, './Data/G_clean_heitor')
readr::write_rds(G_blend, './Data/G_blend_heitor')
readr::write_rds(Aclean, './Data/A_clean_heitor')
readr::write_rds(Aclean, './Data/A_clean_heitor')
readr::write_rds(H, './Data/H_heitor')
readr::write_rds(Hinv.sparse, './Data/Hinv_sparse_heitor')
readr::write_rds(means, './Data/adj_mean_OM_heitor')





