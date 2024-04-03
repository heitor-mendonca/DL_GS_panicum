###############################################################################
#  Article: Genomic selection with allele dosage in Panicum maximum (Jacq.)   #
#           Lara et al., 2018                                                 #
#           Submitted to G3 (Genes|Genomes|Genetics)                          #
#                                                                             #
#  Personal Project : use Deep Learning to perform GS in this article         #
###############################################################################
library(tidyverse)
library(asreml)
library(corrplot)


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

#-----------------------------------------------------
#phenotypic model

model_OM_1 <- asreml(fixed = OM ~ harvest + at(type, 'clone'):genotype +at(type, 'clone'):genotype:harvest , 
                     random = ~ parent + parent:harvest + harvest:block + at(type, 'progeny'):genotype,
                     residual = ~dsum(~units | harvest), #equivalent to idh(harvest):plot
                     workspace = 32e7,
                     data = data)

model_OM_1 <- asreml(fixed = OM ~ harvest + at(type, 'clone'):genotype +at(type, 'clone'):genotype:harvest , 
                     random = ~  parent:harvest + harvest:block + at(type, 'progeny'):genotype:harvest +parent:block:harvest,
                     #residual = ~dsum(~units | harvest), #equivalent to idh(harvest):plot
                     residual = ~id(harvest):plot, 
                     workspace = 32e7,
                     data = data)

wald(model_OM_1)
plot(model_OM_1)
infoCriteria.asreml(model_OM_1) 
summary(model_OM_1, all=T)$varcomp
pred_trait_OM <- predict.asreml(model_OM_1, classify = 'at(type, progeny):genotype', sed=T)
trait_pred_OM <- pred_trait$pvals$predicted.value      