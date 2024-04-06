###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
#                                                                             #
#  Personal Project : use Deep Learning to perform GS in this article         #
###############################################################################

#preparing gentotypic data

library(AGHmatrix)
library(tidyverse)
library(ASRgenomics)
library(asreml)

#load molecular_data.RData as X4

geno <- X4


#getting G 

G <- AGHmatrix::Gmatrix(geno, ploidy = 4,missingValue =NA)


G_sparse <- AGHmatrix::formatmatrix(G, save = F, return = T)



check_G <- kinship.diagnostics(G)
check_G$list.diagonal
check_G$plot.diag

#getting the A matrix

pedigree <- read.table('genealogy.txt', header = T, sep = '\t') #nao usar readr aqui

A <- AGHmatrix::Amatrix(pedigree, ploidy = 4)



#comparing G and A

G2A <- match.G2A(A=A, G=G, clean=TRUE, ord=TRUE, mism=TRUE, RMdiff=TRUE)

G2A$plotG2A

#misaligned



#aligning
Gclean<- G2A$Gclean
Aclean<- G2A$Aclean

G_aligned <- G.tuneup(Gclean,Aclean, align = T) %>% `[[`('Gb')

G_inv <- G.inverse(G_aligned) %>% `[[` ('Ginv')##ill conditioned at first

G_inv_sparse <-AGHmatrix::formatmatrix(G_inv, save = F, return = T) 

G2A <- match.G2A(A=Aclean, G=G_aligned, clean=TRUE, ord=TRUE, mism=TRUE, RMdiff=TRUE)
 

G2A$plotG2A
kinship.diagnostics(G_aligned)


dim(G_aligned)





#obtaining genomic heritability (aligned matrix)

#import means, intersect with G_align


means <- readr::read_delim('adjusted_means_heitor.txt') %>% mutate(genotype = as.character(genotype)) %>% filter(genotype %in% colnames(G_aligned)) %>% 
  mutate(genotype = as.factor(genotype))





str(means)

mod_h2 <- asreml(fixed = predicted.value ~ 1,
                 random = ~ vm(genotype, G_aligned),
                 data = means)


summary(mod_h2)$varcomp
0.1496711 / (0.1496711 + 0.1013527 )


#broad sense h2 0.6930863
#narrow sense h2 0.5962427




#now with the same matrix as lara (blend 99%G 1% A)



G_blend <- G.tuneup(Gclean, Aclean, blend = T, pblend = .01)%>% `[[`('Gb')
dim(G_blend)

G_blend_inv <- G.inverse(G_blend, sparseform = T)$Ginv

mod_h2 <- asreml(fixed = predicted.value ~ 1,
                 random = ~ vm(genotype, G_blend_inv),
                 data = means)


summary(mod_h2)$varcomp

0.1394284/(0.1394284+0.1004284)
#narrow sense h2 0.5962427