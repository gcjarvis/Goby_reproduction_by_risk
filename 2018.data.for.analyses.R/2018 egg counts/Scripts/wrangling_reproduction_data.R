# --------------------------------------
# Description: Getting data set up to send to Dustin
# Author: George C. Jarvis
# Date: Tue Jul 14 11:56:37 2020
# Notes: want to rerun my original analyses for reproduction, but with all the factors, 
#         run them as fixed. Then run model looking at how per capita female body mass may have differed
#         by total density on the reefs. Will try out different things, but likely want to see how female biomass (per capita)
#         is affected by 1) female density, and 2) total density (M + F)
#
#   re: biomass: It seems like we're most interested in female biomass recollected, do fish respond to density by being smaller? What are these environments
# doing to the density and biomass, biomass is interesting. A poor sample, but not a biased sample, is fine. When unbiased data, take everything you can.
#
#   big model, per_capita
#
# 1. start with biomass: what's going on with biomass on each reef 1) in each treatment, and 2) across space and time (trials and years)
#
# want to see if mean mass per reef varied across treatments, density, year, trial (might include final male count and final biomass as well?). Throw everything at it, to see if there were in fact no differences in mass across treatments
# mainly interested in per capita female biomass recollected, but will also also want to see total per capita biomass?
# --------------------------------------

rm(list=ls())

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(emmeans) #for generating least-squares adjusted means from models (will likely help when I have to
## back-transform means if data transformation is needed)
library(car)

library(tidyverse)

#importing dataset####
repro<-read.csv("Data/15.7.20.Jarvis_egg_counts-per_capita_biomass_fem.csv")
repro<-na.omit(repro) #remove cases where no females were recollected

View(repro)
#intial data viz
pairs(repro)

#data manipulation####

#making Year and Trial factors
repro$Year<- as.factor(repro$Year)
repro$Trial<-as.factor(repro$Trial)

#sqrt-transforming egg counts to better satisfy assumptions
repro$sqrt.egg.week<-sqrt(repro$egg.week)

#subsetting data from Trial 6 for comparison of caged vs. uncaged treatments
#subset trial
repro.t6<-repro[repro$Trial==6,]
#subset treatments, caged and uncaged only to test for cage effects
## Note that Treatment factor for trial 6 models should be replaced with "T6.comparison"
repro.t6<-repro.t6[repro.t6$T6.comparison == "High" | repro.t6$T6.comparison == "Uncaged", ]

#exporting wrangled data for those that are not working with these data in R
#Data from Trials 1-5
#write.csv(repro,"Data\\egg_counts_after_data_wrangling.2020.4.7.csv", row.names = FALSE)
#Data from Trial 6 only
#write.csv(repro.t6,"Data\\Cage_effects_egg_counts_after_data_wrangling.2020.4.7.csv", row.names = FALSE)

#15_7_20_analyses####

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#biomass####

head(repro)

#untransformed data
me<-lm(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year*Trial*per_capita_female_biomass*recollection_female, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me)

#sqrt-transformed data look worse than untransformed data?
me<-lm(sqrt(reproduction_per_capita_female_biomass_per_day_per_gram) ~ Treatment*Year*Trial*per_capita_female_biomass*recollection_female, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me)

boxplot(reproduction_per_capita_female_biomass_per_day_per_gram~Trial,repro)

#biomass per reef, including density, total density most likely

boxplot(per_capita_female_biomass~Treatment,repro)
boxplot(per_capita_female_biomass~Treatment*Trial,repro) #seems like I recollected similar biomass within each trial, but there were
# differences among trials
boxplot(per_capita_female_biomass~Trial,repro) #definitely differences by trial, particularly in Trials 4-6

be<-lm(per_capita_female_biomass ~ Treatment*Year*Trial*recollection_female*recollection_male_and_female*recollection_Male, repro)
hist(resid(be))
qqnorm(resid(be))
qqline(resid(be))
plot(be)
summary(be)
anova(be)

be2<-update(be, .~. -(Treatment:Year:recollection_female:recollection_Male))
summary(be2)
anova(be2)

be3<-update(be2, .~. -(Treatment:Trial:recollection_female:recollection_Male))
summary(be3)
anova(be3)

be4<-update(be3, .~. -(Treatment:Year:recollection_male_and_female:recollection_Male))
summary(be4)
anova(be4)

be5<-update(be4, .~. -(Treatment:Trial:recollection_male_and_female:recollection_Male))
summary(be5)
anova(be5)

be6<-update(be5, .~. -(Treatment:recollection_female:recollection_male_and_female:recollection_Male))
summary(be6)
anova(be6)

be7<-update(be6, .~. -(Year:recollection_female:recollection_male_and_female:recollection_Male))
summary(be7)
anova(be7)

be8<-update(be7, .~. -(Trial:recollection_female:recollection_male_and_female:recollection_Male))
summary(be8)
anova(be8)

be9<-update(be8, .~. -(Treatment:Year:Trial:recollection_female:recollection_Male))
summary(be9)
anova(be9)

be10<-update(be9, .~. -(Treatment:Year:Trial:recollection_male_and_female:recollection_Male))
summary(be10)
anova(be10)

be11<-update(be10, .~. -(Treatment:Year:recollection_female:recollection_male_and_female:recollection_Male))
summary(be11)
anova(be11)

be12<-update(be11, .~. -(Year:Trial:recollection_female:recollection_male_and_female:recollection_Male))
summary(be12)
anova(be12)

be13<-update(be12, .~. -(Treatment:Trial:recollection_female:recollection_male_and_female:recollection_Male))
summary(be13)
anova(be13)

be14<-update(be13, .~. -(Treatment:Year:Trial:recollection_female:recollection_male_and_female:recollection_Male))
summary(be14)
anova(be14)

be15<-update(be14, .~. -(Treatment:Trial:recollection_female:recollection_male_and_female))
summary(be15)
anova(be15)

be16<-update(be15, .~. -(Treatment:Year:Trial:recollection_female:recollection_male_and_female))
summary(be16)
anova(be16)

be17<-update(be16, .~. -(Treatment:Year:recollection_female:recollection_male_and_female))
summary(be17)
anova(be17)

be18<-update(be17, .~. -(recollection_female:recollection_male_and_female:recollection_Male))
summary(be18)
anova(be18)

be19<-update(be18, .~. -(Trial:recollection_male_and_female:recollection_Male))
summary(be19)
anova(be19)

be20<-update(be19, .~. -(Year:Trial:recollection_male_and_female:recollection_Male))
summary(be20)
anova(be20)

be21<-update(be20, .~. -(Year:recollection_male_and_female:recollection_Male))
summary(be21)
anova(be21)

be22<-update(be21, .~. -(Treatment:recollection_male_and_female:recollection_Male))
summary(be22)
anova(be22)

be23<-update(be22, .~. -(Trial:recollection_female:recollection_Male))
summary(be23)
anova(be23)

be24<-update(be23, .~. -(Year:Trial:recollection_female:recollection_Male))
summary(be24)
anova(be24) # Year:recollection_female:recollection_Male P = 0.0545

be25<-update(be24, .~. -(Year:recollection_female:recollection_Male))
summary(be25)
anova(be25)

be26<-update(be25, .~. -(Treatment:recollection_female:recollection_Male))
summary(be26)
anova(be26)

be27<-update(be26, .~. -(Trial:recollection_female:recollection_male_and_female))
summary(be27)
anova(be27)

be28<-update(be26, .~. -(Year:Trial:recollection_female:recollection_male_and_female))
summary(be28)
anova(be28)

be29<-update(be28, .~. -(Trial:recollection_female:recollection_male_and_female))
summary(be29)
anova(be29)

be30<-update(be29, .~. -(Year:recollection_female:recollection_male_and_female))
summary(be30)
anova(be30)

be31<-update(be30, .~. -(Treatment:recollection_female:recollection_male_and_female))
summary(be31)
anova(be31)

be32<-update(be31, .~. -(Treatment:Trial:recollection_male_and_female))
summary(be32)
anova(be32)

be33<-update(be32, .~. -(Treatment:Trial:recollection_Male))
summary(be33)
anova(be33)

be34<-update(be33, .~. -(Treatment:Year:Trial:recollection_male_and_female))
summary(be34)
anova(be34)

be35<-update(be34, .~. -(Treatment:Year:Trial:recollection_Male))
summary(be35)
anova(be35)

be36<-update(be35, .~. -(Treatment:Year:recollection_male_and_female))
summary(be36)
anova(be36)

be37<-update(be36, .~. -(Treatment:Year:recollection_Male))
summary(be37)
anova(be37)

be38<-update(be37, .~. -(Treatment:Trial:recollection_female))
summary(be38)
anova(be38)

be39<-update(be38, .~. -(Treatment:Year:Trial:recollection_female))
summary(be39)
anova(be39)

be40<-update(be39, .~. -(Treatment:Year:recollection_female))
summary(be40)
anova(be40)

be41<-update(be40, .~. -(recollection_male_and_female:recollection_Male))
summary(be41)
anova(be41)

be42<-update(be41, .~. -(recollection_female:recollection_Male))
summary(be42)
anova(be42)

be43<-update(be42, .~. -(recollection_female:recollection_male_and_female))
summary(be43)
anova(be43)

be44<-update(be43, .~. -(Trial:recollection_male_and_female))
summary(be44)
anova(be44)

be45a<-update(be44, .~. -(Trial:recollection_Male))
summary(be45a)
anova(be45a)

be45b<-update(be45a, .~. -(Year:Trial:recollection_male_and_female))
summary(be45b)
anova(be45b)

be45c<-update(be45b, .~. -(Year:Trial:recollection_Male))
summary(be45c)
anova(be45c)
#not sure why, but there were some 3-way interactions that snuck their way back into the model? Thought I removed earlier

be46<-update(be45c, .~. -(Year:recollection_male_and_female))
summary(be46)
anova(be46)

be47<-update(be46, .~. -(Year:recollection_Male))
summary(be47)
anova(be47)

be48<-update(be47, .~. -(Treatment:recollection_male_and_female))
summary(be48)
anova(be48)

be49<-update(be48, .~. -(Treatment:recollection_Male))
summary(be49)
anova(be49)

#want to look at the trial* recollection on female biomass = significant at this point, suggesting that female biomass recollected
boxplot(per_capita_female_biomass~Trial*recollection_female,repro)
boxplot(per_capita_female_biomass~Trial,repro)
boxplot(per_capita_female_biomass~Year,repro)
boxplot(per_capita_female_biomass~Year*recollection_female,repro) #not a lot of cases where I recollected a high number of females
boxplot(per_capita_female_biomass~recollection_female,repro) #interesting, more females I recollected, the higher the per capita biomass
boxplot(per_capita_female_biomass~Treatment,repro)
boxplot(per_capita_female_biomass~Treatment*Trial,repro) #seems like patternswere consistent within each trial, in terms of biomass
boxplot(per_capita_female_biomass~Treatment*Year,repro)
boxplot(per_capita_female_biomass~recollection_male_and_female,repro) # per caopita female biomass tended to increase with increased total density on the reefs

boxplot(recollection_female~Treatment*Trial,repro)

boxplot(recollection_female~Reef,repro)

be50<-update(be49, .~. -(Year:recollection_female))
summary(be50)
anova(be50)

be51<-update(be50, .~. -(Treatment:recollection_female))
summary(be51)
anova(be51)

#no treatment x trial interaction
be52<-update(be51, .~. -(Treatment:Trial))
summary(be52)
anova(be52)

be53<-update(be52, .~. -(Treatment:Year:Trial))
summary(be53)
anova(be53)

be54<-update(be53, .~. -(Treatment:Year))
summary(be54)
anova(be54)

# I feel like I should keep total recollection in? P=0.08, will see if the reduced model is better
be55<-update(be54, .~. -(recollection_male_and_female))
summary(be55)
anova(be55)

#now need to remove recollection_male?
be56<-update(be55, .~. -(recollection_Male))
summary(be56)
anova(be56)

#calculating LS means based on this linear model
emmeans(be56, pairwise~Treatment)

be57<-lm(per_capita_female_biomass ~ Treatment+Year+Trial+
           recollection_female + recollection_female*Trial, repro)
anova(be57)

summary(be57)

# reproduction####

# I don't think it's as simple as using the same model that I used for biomass tests.
# I need to include per capita female biomass in the model as well

# did a little bit of model selection, and it seems like per capita reproduction is driven by 
# number of fish recollected (male, female, total), and also by per_capita_female_biomass (duh)

re<-lm(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year*Trial*recollection_female*recollection_male_and_female*recollection_Male*per_capita_female_biomass, repro)
hist(resid(re))
qqnorm(resid(re))
qqline(resid(re))
plot(re)
summary(re)
anova(re)

re2<-update(re, .~. -(Treatment:Year:recollection_female:recollection_Male))
summary(re2)
anova(re2)

re3<-update(re2, .~. -(Treatment:Trial:recollection_female:recollection_Male))
summary(re3)
anova(re3)

re4<-update(re3, .~. -(Treatment:Year:recollection_male_and_female:recollection_Male))
summary(re4)
anova(re4)

re5<-update(re4, .~. -(Treatment:Trial:recollection_male_and_female:recollection_Male))
summary(re5)
anova(re5)

re6<-update(re5, .~. -(Treatment:recollection_female:recollection_male_and_female:recollection_Male))
summary(re6)
anova(re6)

re7<-update(re6, .~. -(Year:recollection_female:recollection_male_and_female:recollection_Male))
summary(re7)
anova(re7)

re8<-update(re7, .~. -(Trial:recollection_female:recollection_male_and_female:recollection_Male))
summary(re8)
anova(re8)

re9<-update(re8, .~. -(Treatment:Year:Trial:recollection_female:recollection_Male))
summary(re9)
anova(re9)

re10<-update(re9, .~. -(Treatment:Year:Trial:recollection_male_and_female:recollection_Male))
summary(re10)
anova(re10)

re11<-update(re10, .~. -(Treatment:Year:recollection_female:recollection_male_and_female:recollection_Male))
summary(re11)
anova(re11)

re12<-update(re11, .~. -(Year:Trial:recollection_female:recollection_male_and_female:recollection_Male))
summary(re12)
anova(re12)

re13<-update(re12, .~. -(Treatment:Trial:recollection_female:recollection_male_and_female:recollection_Male))
summary(re13)
anova(re13)

re14<-update(re13, .~. -(Treatment:Year:Trial:recollection_female:recollection_male_and_female:recollection_Male))
summary(re14)
anova(re14)

re15<-update(re14, .~. -(Treatment:Trial:recollection_female:recollection_male_and_female))
summary(re15)
anova(re15)

re16<-update(re15, .~. -(Treatment:Year:Trial:recollection_female:recollection_male_and_female))
summary(re16)
anova(re16)

re17<-update(re16, .~. -(Treatment:Year:recollection_female:recollection_male_and_female))
summary(re17)
anova(re17)

re18<-update(re17, .~. -(recollection_female:recollection_male_and_female:recollection_Male))
summary(re18)
anova(re18)

re19<-update(re18, .~. -(Trial:recollection_male_and_female:recollection_Male))
summary(re19)
anova(re19)

re20<-update(re19, .~. -(Year:Trial:recollection_male_and_female:recollection_Male))
summary(re20)
anova(re20)

re21<-update(re20, .~. -(Year:recollection_male_and_female:recollection_Male))
summary(re21)
anova(re21)

re22<-update(re21, .~. -(Treatment:recollection_male_and_female:recollection_Male))
summary(re22)
anova(re22)

re23<-update(re22, .~. -(Trial:recollection_female:recollection_Male))
summary(re23)
anova(re23)

re24<-update(re23, .~. -(Year:Trial:recollection_female:recollection_Male))
summary(re24)
anova(re24) # Year:recollection_female:recollection_Male P = 0.0545

re25<-update(re24, .~. -(Year:recollection_female:recollection_Male))
summary(re25)
anova(re25)

re26<-update(re25, .~. -(Treatment:recollection_female:recollection_Male))
summary(re26)
anova(re26)

re27<-update(re26, .~. -(Trial:recollection_female:recollection_male_and_female))
summary(re27)
anova(re27)

re28<-update(re26, .~. -(Year:Trial:recollection_female:recollection_male_and_female))
summary(re28)
anova(re28)

re29<-update(re28, .~. -(Trial:recollection_female:recollection_male_and_female))
summary(re29)
anova(re29)

boxplot(per_capita_female_biomass~Trial*recollection_female,repro)

re30<-update(re29, .~. -(Year:recollection_female:recollection_male_and_female))
summary(re30)
anova(re30)

re31<-update(re30, .~. -(Treatment:recollection_female:recollection_male_and_female))
summary(re31)
anova(re31)

re32<-update(re31, .~. -(Treatment:Trial:recollection_male_and_female))
summary(re32)
anova(re32)

re33<-update(re32, .~. -(Treatment:Trial:recollection_Male))
summary(re33)
anova(re33)

re34<-update(re33, .~. -(Treatment:Year:Trial:recollection_male_and_female))
summary(re34)
anova(re34)

re35<-update(re34, .~. -(Treatment:Year:Trial:recollection_Male))
summary(re35)
anova(re35)

re36<-update(re35, .~. -(Treatment:Year:recollection_male_and_female))
summary(re36)
anova(re36)

re37<-update(re36, .~. -(Treatment:Year:recollection_Male))
summary(re37)
anova(re37)

re38<-update(re37, .~. -(Treatment:Trial:recollection_female))
summary(re38)
anova(re38)

re39<-update(re38, .~. -(Treatment:Year:Trial:recollection_female))
summary(re39)
anova(re39)

re40<-update(re39, .~. -(Treatment:Year:recollection_female))
summary(re40)
anova(re40)

re41<-update(re40, .~. -(recollection_male_and_female:recollection_Male))
summary(re41)
anova(re41)

re42<-update(re41, .~. -(recollection_female:recollection_Male))
summary(re42)
anova(re42)

re43<-update(re42, .~. -(recollection_female:recollection_male_and_female))
summary(re43)
anova(re43)

re44<-update(re43, .~. -(Trial:recollection_male_and_female))
summary(re44)
anova(re44)

re45a<-update(re44, .~. -(Trial:recollection_Male))
summary(re45a)
anova(re45a)

re45b<-update(re45a, .~. -(Year:Trial:recollection_male_and_female))
summary(re45b)
anova(re45b)

re45c<-update(re45b, .~. -(Year:Trial:recollection_Male))
summary(re45c)
anova(re45c)
#not sure why, but there were some 3-way interactions that snuck their way back into the model? Thought I removed earlier

re46<-update(re45c, .~. -(Year:recollection_male_and_female))
summary(re46)
anova(re46)

re47<-update(re46, .~. -(Year:recollection_Male))
summary(re47)
anova(re47)

re48<-update(re47, .~. -(Treatment:recollection_male_and_female))
summary(re48)
anova(re48)

re49<-update(re48, .~. -(Treatment:recollection_Male))
summary(re49)
anova(re49)

re50<-update(re49, .~. -(Year:recollection_female))
summary(re50)
anova(re50)

re51<-update(re50, .~. -(Treatment:recollection_female))
summary(re51)
anova(re51)

#no treatment x trial interaction
re52<-update(re51, .~. -(Treatment:Trial))
summary(re52)
anova(re52)

re53<-update(re52, .~. -(Treatment:Year:Trial))
summary(re53)
anova(re53)

re54<-update(re53, .~. -(Treatment:Year))
summary(re54)
anova(re54)

# I feel like I should keep total recollection in? P=0.08, will see if the reduced model is retter
re55<-update(re54, .~. -(recollection_male_and_female))
summary(re55)
anova(re55)

#now need to remove recollection_male?
re56<-update(re55, .~. -(recollection_Male))
summary(re56)
anova(re56)

#calculating LS means based on this linear model
emmeans(re56, pairwise~Treatment)

# trying to run nested model with fixed factors#####
summary(repro)
head(repro)

nest1<-aov(per_capita_female_biomass~Treatment%in%Trial%in%Year+Year+((Trial%in%Year)*recollection_female),data = repro) 
summary(nest1)
anova(nest1)
