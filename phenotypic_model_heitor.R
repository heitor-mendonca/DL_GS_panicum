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


# Checking trait data

nutritional_traits <- readr::read_delim('Nutritional_traits.txt')


nutritional_traits %>% select(OM) %>% drop_na() %>% boxplot()
nutritional_traits %>% select(CP) %>% drop_na() %>% boxplot()
nutritional_traits %>% select(IVD) %>% drop_na() %>% boxplot()


productive_traits <- readr::read_delim('Productive_traits.txt')


productive_traits %>% select(LDM) %>% drop_na() %>% boxplot()

productive_traits %>% select(PLB) %>% drop_na() %>% boxplot()


plot(density(productive_traits$PLB%>%na.omit()))
plot(density(productive_traits$LDM%>%na.omit()))
plot(density(nutritional_traits$OM%>%na.omit()))
plot(density(nutritional_traits$CP%>%na.omit()))
plot(density(nutritional_traits$IVD%>%na.omit()))

#OM seems to be the best, we're going to work with it for the moment
#---------------------------------------------------

data <- readr::read_delim('Nutritional_traits.txt') %>% select(!c(CP,IVD)) %>%mutate(genotype = as.factor(genotype),
                                                                                     type     = as.factor(type),
                                                                                     parent   = as.factor(parent),
                                                                                     block    = as.factor(block),
                                                                                     harvest  = as.factor(harvest),
                                                                                     plot     = as.factor(plot))
str(data)
head(data)

#----------------------------------------------------
#exploratory analyses


ggplot(data, aes(y = OM)) + 
  geom_boxplot()+
  theme_minimal()


ggplot(data, aes(x = OM) + 
  geom_density()+
  theme_minimal()

summary(data)


data %>% drop_na() %>% summarise(mean = mean(OM), .by = harvest) %>% ggplot(aes(x = harvest, y = mean))+
  geom_point()

data %>% drop_na() %>% group_by(genotype, harvest) %>% summarise(mean = mean(OM)) %>% ggplot(aes(x = harvest, y = mean, group = genotype)) + geom_line()


data %>% fct_lump_prop(genotypes,prop = 4) %>% drop_na()  %>% count(genotype) %>% arrange(n)
#-----------------------------------------------------
#phenotypic model

model_OM <- asreml(fixed = OM ~ harvest + at(type,'clone'):genotype +at(type, 'clone'):genotype:harvest  , 
                     random = ~ parent:block:harvest  + harvest:block +  parent:harvest+ 
                     at(type, 'progeny'):(genotype):exp(harvest),
                     residual = ~corh(harvest):plot, 
                     workspace = 32e7,
                     data = data)


model_OM <- update(model_OM)

wald(model_OM)
plot(model_OM)
infoCriteria.asreml(model_OM) 
summary(model_OM, all=T)$varcomp


#---------------------------------------------------------------------------------------------------------------------------------
model_OM_h2 <- asreml(fixed = OM ~ harvest  , 
                   random = ~ parent:block:harvest  + harvest:block +  parent:harvest+ 
                     genotype:exp(harvest),
                   residual = ~corh(harvest):plot, 
                   workspace = 32e7,
                   data = data)



summary(model_OM_h2, all=T)$varcomp

vc.g <- 0.50472274
vc.g

# Mean variance of a difference of two genotypic BLUPs
vdBLUP.mat <- predict.asreml(model_OM_h2, classify="genotype", sed=TRUE)$sed^2 # obtain squared s.e.d. matrix 
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
vdBLUP.avg #0.05455038
H2Cullis <- 1 - (vdBLUP.avg / 2 / vc.g)
H2Cullis #0.8091336
#--------------------------------------------------------------------------------------------



pred_trait_OM <- predict.asreml(model_OM_1, classify = 'at(type, progeny):genotype', sed=T, pworkspace = 32e7)
pred_trait_OM <- predict.asreml(model_OM_1, classify = 'genotype', pworkspace = 32e7)
trait_pred_OM <- pred_trait_OM$pvals$predicted.value      
