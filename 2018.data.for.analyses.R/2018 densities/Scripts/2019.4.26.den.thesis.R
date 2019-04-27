# Description: Density analyses for thesis, decided to go with day, not week, so let's see if this has a big effect
# Author: George C Jarvis
# Date: Fri Apr 26 21:39:20 2019
# --------------

#I think if I go with day for all analyses, it's better

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(MASS)
library(nlme)
library(pwr)

gob.den<-read.csv("Data/density.2019.4.26.csv")
gob.den<- na.omit(gob.den) #remove NA's

#subsetting by trial
egg.2017<-egob.den[(gob.den$Trial>3) & (gob.den$Trial<6), ]
egg.2018.t4.5<-gob.den[(gob.den$Trial>3) & (gob.den$Trial<6), ]
egg.2018.t4.5<-gob.den[(gob.den$Trial>3) & (gob.den$Trial<6), ]
