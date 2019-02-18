
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
                    genotype=factor(genotype), plot=factor(plot))
str(traits)

trait <- asreml(fixed = OM ~ harvest + at(type,1):genotype + at(type,1):genotype:harvest,
                random = ~ block:harvest + parent:block:harvest +
                #          at(type,2):genotype:id(harvest),
                #          at(type,2):genotype:diag(harvest),
                #          at(type,2):genotype:cor(harvest),
                #          at(type,2):genotype:corh(harvest),
                #          at(type,2):genotype:ar1(harvest),
                #          at(type,2):genotype:ar1h(harvest),
                #          at(type,2):genotype:exp(harvest),
                          at(type,2):genotype:exph(harvest),
                #          at(type,2):genotype:us(harvest),
                #rcov = ~ id(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ diag(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ cor(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ corh(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ ar1(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                rcov = ~ ar1h(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ exp(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ exph(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ us(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                data = traits) 
wald(trait)
plot(trait)
info.crit.asreml(trait) 
summary(trait, all=T)$varcomp
pred_trait <- predict(trait, classify = 'at(type, progeny):genotype', sed=T)
trait_pred <- pred_trait$predictions$pvals$predicted.value            

# Means
offspring <- pred_trait$predictions$pvals$genotype[574:1143]
trait_means <- cbind(offspring, trait_pred[574:1143])
head(trait_means); tail(trait_means)
trait_means <- round(trait_means, digits = 4)

# Generalized heritability
trait.h2 <- asreml(fixed = OM ~ harvest + at(type,1):genotype + at(type,1):genotype:harvest,
                   random = ~ block:harvest + parent:block:harvest + at(type,2):genotype + at(type,2):genotype:harvest,
                   rcov = ~ id(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                   data = traits) 
info.crit.asreml(trait.h2) 
summary(trait.h2, all=T)$varcomp
pred2_trait.h2 <- predict(trait.h2, classify = 'at(type, progeny):genotype')
pev_trait <- (pred2_trait.h2$predictions$avsed)^2                              
h2.trait <- 1-(pev_trait/(2*summary(pred2_trait.h2, all=T)$varcomp[3,2]))      


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
                random = ~ block:harvest + parent:block:harvest +
                #          at(type,2):genotype:id(harvest),
                #          at(type,2):genotype:diag(harvest),
                #          at(type,2):genotype:cor(harvest),
                #          at(type,2):genotype:corh(harvest),
                #          at(type,2):genotype:ar1(harvest),
                #          at(type,2):genotype:ar1h(harvest),
                #          at(type,2):genotype:exp(harvest),
                          at(type,2):genotype:exph(harvest),
                #          at(type,2):genotype:us(harvest),
                #rcov = ~ id(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ diag(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ cor(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ corh(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ ar1(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ ar1h(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ exp(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                rcov = ~ exph(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                #rcov = ~ us(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                data = traits) 
wald(trait)
plot(trait)
info.crit.asreml(trait) 
summary(trait, all=T)$varcomp
pred_trait <- predict(trait, classify = 'at(type, progeny):genotype', sed=T)
trait_pred <- pred_trait$predictions$pvals$predicted.value            

# Means
offspring <- pred_trait$predictions$pvals$genotype[574:1143]
trait_means <- cbind(offspring, trait_pred[574:1143])
head(trait_means); tail(trait_means)
trait_means <- round(trait_means, digits = 4)

# Generalized heritability
trait.h2 <- asreml(fixed = LDM ~ harvest + at(type,1):genotype + at(type,1):genotype:harvest,
                   random = ~ block:harvest + parent:block:harvest + at(type,2):genotype + at(type,2):genotype:harvest,
                   rcov = ~ id(harvest):plot, maxiter=900, control = asreml.control(workspace=32e7), 
                   data = traits) 
info.crit.asreml(trait.h2) 
summary(trait.h2, all=T)$varcomp
pred2_trait.h2 <- predict(trait.h2, classify = 'at(type, progeny):genotype')
pev_trait <- (pred2_trait.h2$predictions$avsed)^2                              
h2.trait <- 1-(pev_trait/(2*summary(pred2_trait.h2, all=T)$varcomp[2,2]))      


######################################################################################################################

