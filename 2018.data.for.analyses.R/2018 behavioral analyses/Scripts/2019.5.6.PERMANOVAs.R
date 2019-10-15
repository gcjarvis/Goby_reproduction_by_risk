# Description: PERMANOVA for behavioral observations
# Author: George C Jarvis
# Date: Mon May 06 10:19:10 2019
# --------------

#data for behaviors likely do not satisfy the pparametric assumptions of normality and
# equal variance, so went with premanovato test for differnces in means or dispersions

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plyr)
library(vegan)
library(psych)
library(vegan)

#subsetting df that includes all trials
behave<-read.csv("Data/2019.4.25.behavior.ss.csv") #doesn't include unnecessary NA's
behave<- na.omit(behave)
behave<-arrange(behave, Treatment) #ordering by treatment for permanova and permdisp
colSums(is.na(behave))# no NA's

#subsetting by year
behave.2017<-behave[behave$Trial<4, ]
View(behave.2017)
behave.2018.t4.5<-behave[c((behave$Trial>3) & (behave$Trial<6)), ]
behave.2018.t4.5<- na.omit(behave.2018.t4.5)

#1)permanova for behaviors in 2017
#setting up object with all response variables
#made separate excel file bc of issues
behave<-read.csv("Data/2019.5.6..behavior.ss.trial.1.3.csv") #doesn't include unnecessary NA's
behave<- na.omit(behave)
behave<-arrange(behave, Treatment) #ordering by treatment for permanova and permdisp
colSums(is.na(behave))# no NA's


object<-cbind(behave.2017$proportion.exposed,behave.2017$movements.min,
              behave.2017$bites.min, behave.2017$total.dist.moved,
              behave.2017$courtship.min)

object<-behave[,4:8]

View(object)

permanovamodel1<-adonis(object~Treatment*Trial, 
                        data = behave, permutations = 100,
                        method="euclidean", by= "terms")
permanovamodel1 #suggests no differences in multivariate space, and no differences
#in behaviors by trial

#permdisp for counts
## euclidean distances between samples
dis1 <- vegdist(object, method = "euclidean")

#including each factor, along with the number of replicates for each (#factor, #replicates)

groups <- factor(c(rep(1,18),rep(2,18),rep(3,18)), labels =
                   c("High","Low","Medium"))

groups
dis1

## Calculate multivariate dispersions
mod1 <- betadisper(dis1, groups)
mod1

## Perform test
anova(mod1)

## Permutation test for F
permutest(mod1, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod1))
plot(mod.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod1)

#2017, doesn't appear to be any differences in behavior
#among treatments, no differences in location of means (permanova results insig.)
#and no differences in the dispersion of points about each mean (permdisp insig.)

#2)permanova for behaviors in 2018
#setting up object with all response variables
#made separate excel file bc of issues
behave2<-read.csv("Data/2019.5.6..behavior.ss.trial.4.5.csv") #doesn't include unnecessary NA's
behave2<- na.omit(behave2)
behave2<-arrange(behave2, Treatment) #ordering by treatment for permanova and permdisp
colSums(is.na(behave2))# no NA's
View(behave2)

#object<-cbind(behave.2017$proportion.exposed,behave.2017$movements.min,
#              behave.2017$bites.min, behave.2017$total.dist.moved,
#              behave.2017$courtship.min)

object2<-behave2[,4:8]

View(object2)

permanovamodel2<-adonis(object2~Treatment*Trial, 
                        data = behave2, permutations = 100,
                        method="euclidean", by= "terms")
permanovamodel2 #suggests no differences in multivariate space, and no differences
#in behaviors by trial

#permdisp for counts
## euclidean distances between samples
dis2 <- vegdist(object2, method = "euclidean")

#including each factor, along with the number of replicates for each (#factor, #replicates)
groups2 <- factor(c(rep(1,12),rep(2,12),rep(3,11)), labels =
                   c("High","Low","Medium"))

View(behave2)
#groups
#dis2

## Calculate multivariate dispersions
mod2 <- betadisper(dis2, groups2)
mod2

## Perform test
anova(mod2)

## Permutation test for F
permutest(mod2, pairwise = TRUE)

## Tukey's Honest Significant Differences
mod.HSD2 <- TukeyHSD(mod2)
plot(mod.HSD2)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod2)

#it seems like there was an effect of trial on the behaviors
#observed, so I'm going redo with df ordered by trial

behave3<-read.csv("Data/2019.5.6..behavior.ss.trial.4.5.csv") #doesn't include unnecessary NA's
behave3<- na.omit(behave3)
behave3<-arrange(behave3, Trial) #ordering by treatment for permanova and permdisp
colSums(is.na(behave3))# no NA's

View(behave3)
object3<-behave3[,4:8]

View(object3)

permanovamodel1<-adonis(object3~Treatment*Trial, 
                        data = behave3, permutations = 100,
                        method="euclidean", by= "terms")
permanovamodel1 #suggests no differences in multivariate space, and no differences
#in behaviors by trial

#permdisp for counts
## euclidean distances between samples
dis3 <- vegdist(object3, method = "euclidean")

#including each factor, along with the number of replicates for each (#factor, #replicates)
groups3 <- factor(c(rep(1,17),rep(2,18)), labels =
                   c("4","5"))

groups3
dis3

## Calculate multivariate dispersions
mod3 <- betadisper(dis3, groups3)
mod3

## Perform test
anova(mod3)

## Permutation test for F
permutest(mod3, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod3))
plot(mod.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod3)

