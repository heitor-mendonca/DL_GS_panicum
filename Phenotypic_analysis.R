
###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
###############################################################################


rm(list=ls());ls()
library(asreml)
library(asremlPlus)
library(dae)
library(dplyr)

####################
# Nutritional values
traits <- read.table('Nutritional_traits.txt', h=T, sep='\t', dec='.')
str(traits)
head(traits, 8)

traits <- transform(traits, harvest=factor(harvest), block=factor(block), parent=factor(parent), 
                    genotype=factor(genotype), plot=factor(plot), type = factor(type))
str(traits)

trait <- asreml(fixed = OM ~ harvest + at(type,1):genotype + at(type,1):genotype:harvest,
                random = ~ block:harvest + parent:harvest + parent:block:harvest +
                          at(type,2):genotype:id(harvest),
                #          at(type,2):genotype:diag(harvest),
                #          at(type,2):genotype:cor(harvest),
                #          at(type,2):genotype:corh(harvest),
                #          at(type,2):genotype:ar1(harvest),
                #          at(type,2):genotype:ar1h(harvest),
                #          at(type,2):genotype:exp(harvest),
                #          at(type,2):genotype:exph(harvest),
                #          at(type,2):genotype:us(harvest),
                residual = ~ id(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ diag(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ cor(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ corh(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ ar1(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ ar1h(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ exp(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ exph(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ us(harvest):plot, maxiter=900, workspace=32e7, 
                data = traits) 
wald(trait)
plot(trait)
infoCriteria.asreml(trait) 
summary(trait, all=T)$varcomp
pred_trait <- predict.asreml(trait, classify = 'at(type, progeny):genotype', sed=T)
trait_pred <- pred_trait$pvals$predicted.value            

# Means
offspring <- pred_trait$pvals$genotype[574:1143]
trait_means <- cbind(offspring, trait_pred[574:1143])
head(trait_means); tail(trait_means)
trait_means <- round(trait_means, digits = 4)

# Generalized heritability
trait.h2 <- asreml(fixed = OM ~ harvest + at(type,1):genotype + at(type,1):genotype:harvest,
                   random = ~ block:harvest + parent:harvest + parent:block:harvest + at(type,2):genotype + at(type,2):genotype:harvest,
                   residual = ~ id(harvest):plot, maxiter=900, workspace=32e7, 
                   data = traits) 
infoCriteria.asreml(trait.h2) 
summary(trait.h2, all=T)$varcomp
pred2_trait.h2 <- predict.asreml(trait.h2, classify = 'at(type, progeny):genotype')
pev_trait <- (pred2_trait.h2$avsed)^2                              
h2.trait <- 1-(pev_trait/(2*summary(trait.h2, all=T)$varcomp[4,1]))      


######################################################################################################################

rm(list=ls());ls()
####################
# Nutritional values
traits <- read.table('Productive_traits.txt', h=T, sep='\t', dec='.')
str(traits)
head(traits, 8)

traits <- transform(traits, harvest=factor(harvest), block=factor(block), parent=factor(parent), 
                    genotype=factor(genotype), plot=factor(plot))
str(traits)

trait <- asreml(fixed = LDM ~ harvest + at(type,1):genotype + at(type,1):genotype:harvest,
                random = ~ block:harvest + parent:harvest + parent:block:harvest +
                          at(type,2):genotype:id(harvest),
                #          at(type,2):genotype:diag(harvest),
                #          at(type,2):genotype:cor(harvest),
                #          at(type,2):genotype:corh(harvest),
                #          at(type,2):genotype:ar1(harvest),
                #          at(type,2):genotype:ar1h(harvest),
                #          at(type,2):genotype:exp(harvest),
                #          at(type,2):genotype:exph(harvest),
                #          at(type,2):genotype:us(harvest),
                residual = ~ id(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ diag(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ cor(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ corh(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ ar1(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ ar1h(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ exp(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ exph(harvest):plot, maxiter=900, workspace=32e7, 
                #residual = ~ us(harvest):plot, maxiter=900, workspace=32e7, 
                data = traits) 
wald(trait)
plot(trait)
infoCriteria.asreml(trait) 
summary(trait, all=T)$varcomp
pred_trait <- predict.asreml(trait, classify = 'at(type, progeny):genotype', sed=T)
trait_pred <- pred_trait$pvals$predicted.value            

# Means
offspring <- pred_trait$pvals$genotype[574:1143]
trait_means <- cbind(offspring, trait_pred[574:1143])
head(trait_means); tail(trait_means)
trait_means <- round(trait_means, digits = 4)

# Generalized heritability
trait.h2 <- asreml(fixed = LDM ~ harvest + at(type,1):genotype + at(type,1):genotype:harvest,
                   random = ~ block:harvest + parent:harvest + parent:block:harvest + at(type,2):genotype + at(type,2):genotype:harvest,
                   residual = ~ id(harvest):plot, maxiter=900, workspace=32e7, 
                   data = traits) 
infoCriteria.asreml(trait.h2) 
summary(trait.h2, all=T)$varcomp
pred2_trait.h2 <- predict.asreml(trait.h2, classify = 'at(type, progeny):genotype')
pev_trait <- (pred2_trait.h2$avsed)^2                              
h2.trait <- 1-(pev_trait/(2*summary(trait.h2, all=T)$varcomp[3,1]))      


######################################################################################################################

