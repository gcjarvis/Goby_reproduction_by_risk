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
ptl<-read.csv("Data/2019.3.19.ptl.csv")
#arranging by treatment for PERMANOVA + PermDisp
ptl<-arrange(ptl, Treatment)

#subsetting data 
d2<-subset(ptl,Trial<6) #trial 4 and 5
count(d2, c("Reef","Treatment"))#what is the sample size for each treatment
d3<-subset(ptl,Trial>5) #trial 6

#need to get rid of NA's to run PERMANOVA
colSums(is.na(d2)) # count of # of NA's within df
#omit na within score column
d2<-d2[!is.na(d2$Score.zero),]
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

permanovamodel1<-adonis(Score.zero~Treatment, data = d3, permutations = 100, method="euclidean", by= "terms")
permanovamodel1
