# Description: testing out MCMC glmm package with model from before
# Author: George C Jarvis
# Date: Sun Mar 08 16:55:19 2020
# Notes: not sure if this will work, took code and notes from 
# https://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
# --------------

#load package
library()

#loading data####
repro<-read.csv("Data/new.data.2019.9.30.csv")
repro<-na.omit(repro) # no NA's to omit

#data manipulation####
#adding column for average density, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
repro$Year <- ifelse(repro$Trial <=3, 2017, 2018)
#want it as a factor? Going to make a variable with it as a factor, run the model, and see if I get different results
repro$Year.fact<- as.factor(repro$Year)

#adding egg/week variable to the dataset
repro<-repro %>%
  mutate(egg.week = ifelse(Trial<4, Egg.count/1,
                           ifelse(Trial == 4| Trial == 5, (ceiling(Egg.count/4)),
                                  ifelse(Trial == 6, (ceiling(Egg.count/2)), NA))))

#modeling

#pick up here next time after emailing Steve to see what he has to say

library(MCMCglmm)
## Loading required package: coda
## Loading required package: lattice
## Loading required package: ape

#no idea if these priors are correct for my model
prior = list(R = list(V = 1, n = 0, fix = 1), G = list(G1 = list(V = 1, n = 1),
         G2 = list(V = 1, n = 1), G3 = list(V = 1, n = 1), G4 = list(V = 1, n = 1),
         G5 = list(V = 1, n = 1)))
set.seed(45)
MCMC <- MCMCglmm(profit ~ 1, random = ~year + farmer + place + gen + district,
                 data = farmers, family = "categorical", prior = prior, verbose = FALSE)
summary(MCMC)
