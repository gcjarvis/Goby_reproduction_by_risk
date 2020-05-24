# Description: Script for recollctions (i.e. survival) from Jarvis and Steele
# Author: George C Jarvis
# Date: Sun May 24 13:48:41 2020
# Notes:
# --------------

rm(list=ls())

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(emmeans) #for generating least-squares adjusted means from models

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

# import dataset ####
reco<-read.csv("Data/2019.10.8.recollection.data.csv")

# data manipulation ####

#adding column for survivorship, dividing recollections by 20 (initial number of fish on reefs)
reco$Survivorship<-reco$Count/20

#changing trial and year to factors
reco$Trial<- as.factor(reco$Trial)
reco$Year.fact<- as.factor(reco$Year)


#subsetting data from trial 6 to compare recollections between HR caged and uncaged treatments
reco.t6<-reco[reco$Trial==6,]

#subset treatments, caged and uncaged only (T6.comparison)
#subsetting by multiple character factors within a variable. Awesome code!
reco.t6.comp<-reco.t6[reco.t6$T6.comparison == "High" | reco.t6$T6.comparison == "Control", ]

#only need to do a t-test for trial 6
# are there differences in survivorship between caged and uncaged reefs?

# analyses ####

# analyzing mixed models with log-likelihood estimates and chi-square tests.
# The full model for trials 1-5 will include treatment, year, the interaction (all fixed), 
## plus trial nested within year, and treatment x trial nested within year as random effects.
# I am planning to keep Trial in the model, but treatment x trial will likely be removed if it is not significant

# Start with full model, then reduce model first by non-significant random effects, then by NS. fixed effects
# want to definitely leave in Trial term as random effect, and Treatement, Year, T x Y, and avg.inhab as fixed effects.
# remove all NS. interactions with the covariate (fixed and random)



