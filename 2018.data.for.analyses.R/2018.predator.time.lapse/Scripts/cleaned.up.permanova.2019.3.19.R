# Description: Permanova for predator time lapses
# Author: George C Jarvis
# Date: Tue Mar 19 18:34:17 2019
# --------------

#NOTE: these data only incorporate PTL reads up through 2019.3.19

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

#importing raw data
#old df ptl<-read.csv("Data/2019.3.19.ptl.csv")
#new df as of 2019.4.26
ptl<-read.csv("Data/2019.4.26.ptl.permanova.csv")
#arranging by treatment for PERMANOVA + PermDisp
ptl<-arrange(ptl, Treatment)
View(ptl)
colSums(is.na(ptl))# no NA's

#after lunch: get results from permanova and make tables (I'd make a new script for this
# copying what you did from this line and above)

#after that, make rest of figures and tables in ggplot and 

#after that, finish draft...


#subsetting data 
d2<-subset(ptl,Trial<6) #trial 4 and 5
count(d2, c("Reef","Treatment"))#what is the sample size for each treatment
d3<-subset(ptl,Trial>5) #trial 6

#need to get rid of NA's to run PERMANOVA
colSums(is.na(d2)) # count of # of NA's within df
#omit na within score column
d2<-d2[!is.na,]
#removing control treatment from df for t4 and 5 (will need to do this for 2017 data as well)
d2<-subset(ptl,Treatment!="Control")

toBeRemoved<-which(d2$Treatment=="Control")
d2<-d2[-toBeRemoved,]

#same thing for d3
colSums(is.na(d3)) # count of # of NA's within df
#omit na within score column
d3<-d3[!is.na(d3$Score.zero),]

#PERMANOVA, one-way by Treatment, with score and count as two separate models
#using data with zeros
#t4 and 5
#Score

#figuring out how to get rid of the control

permanovamodel1<-adonis(contain.predator~Treatment, 
                        data = ptl, permutations = 100,
                        method="euclidean", by= "terms")
permanovamodel1

#permdisp
## Bray-Curtis distances between samples
dis <- vegdist(ptl$contain.predator)

#dis <- vegdist(varespec)

## First 16 sites grazed, remaining 8 sites ungrazed
groups <- factor(c(rep(1,16), rep(2,8)), labels = c("grazed","ungrazed"))

## Calculate multivariate dispersions
mod <- betadisper(dis, groups)
mod

## Perform test
anova(mod)

## Permutation test for F
permutest(mod, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)

## simulate missing values in 'd' and 'group'
groups[c(2,20)] <- NA
dis[c(2, 20)] <- NA
mod2 <- betadisper(dis, groups) ## warnings
mod2
permutest(mod, control = permControl(nperm = 100))
anova(mod2)
plot(mod2)
boxplot(mod2)
plot(TukeyHSD(mod2))