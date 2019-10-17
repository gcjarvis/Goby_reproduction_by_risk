# Description: MANCOVA with data from all trials for manuscript
# Author: George C Jarvis
# Date: Tue Oct 08 16:21:11 2019
# Notes:
# --------------

#clear workspace
rm(list=ls())

library(sciplot)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(plyr)
#packages needed for MANCOVA
library(car)#type III SS
library(psych)#descriptive statistics
library(effects)#adjusted means
#easiest mancova package
library(jmv)

behave<-read.csv("Data/2019.4.25.behavior.ss.csv") #doesn't include unnecessary NA's
head(behave)
#subsetting by year
behave.2017<-behave[behave$Trial<4, ]
behave.2018.t4.5<-behave[c((behave$Trial>3) & (behave$Trial<6)), ]

# FYI, ran both of these MANCOVAs with Trial included as a factor, and didn't find any differences by trial
# removed trial from the model

#MANCOVAs with density#####

#2017
mancova(data=behave.2017, deps=c('proportion.exposed','movements.min','bites.min',
                                 'total.dist.moved','courtship.min'),factors='Treatment',covs = 'density')
#only significant factor was coutrship rates by density
# it seems like there was a positive relationship between the two, but not enough to drie differences
plot(courtship.min~density,data=behave.2017)

#2018
mancova(data=behave.2018.t4.5, deps=c('proportion.exposed','movements.min','bites.min',
                                      'total.dist.moved','courtship.min'),factors='Treatment',covs = 'density')
#phew, saw the same trends as before:
# overall differences, driven by differences in exposure time

#low > med > high: trendlikely driven by differences between low and high treatments (i.e. low and med were very similar)
boxplot(proportion.exposed~Treatment, data=behave.2018.t4.5)

#calculatign percent differences
tapply(behave.2018.t4.5$proportion.exposed, 
       behave.2018.t4.5$Treatment, mean) #18% (sig) decrease in prop time exposed on HR
tapply(behave.2017$proportion.exposed, 
       behave.2017$Treatment, mean) #22% (NS) decrease in prop time exposed on HR


#also saw effect of covariate on 1) prop.exposed and 2) total.dist.moved

#plotting regressions for general idea
# 1) prop exposed by density
plot(proportion.exposed~density,data=behave.2018.t4.5)
# exposure time seems to increase with increasing density, because:
# a) densities were likely higher on low risk treatments --> issue with confounding treatment and density

boxplot(density~Treatment,data=behave.2018.t4.5)#data shows this, though it may not be statistically sig.

#MANCOVAs with recollection as covariate########

#2017
mancova(data=behave.2017, deps=c('proportion.exposed','movements.min','bites.min',
                                 'total.dist.moved','courtship.min'),factors='Treatment',covs = 'recollection')
#no overall, or covariate effects

#2018
mancova(data=behave.2018.t4.5, deps=c('proportion.exposed','movements.min','bites.min',
                                      'total.dist.moved','courtship.min'),factors='Treatment',covs = 'recollection')
#phew, saw the same trends as before:
# overall differences, driven by differences in exposure time, but no effect of recollection on bahaviors overall

#low > med > high: trendlikely driven by differences between low and high treatments (i.e. low and med were very similar)
boxplot(proportion.exposed~Treatment, data=behave.2018.t4.5)

#plotting regressions for general idea
# 1) prop exposed by density
plot(proportion.exposed~recollection,data=behave.2018.t4.5)
# no effect of recollection on behaviors (no real relationship between exposure time and recollections)

# MANOVAs #####
# I think this is the way to go, becuase densities are confounded by treatment, and recollections didn't seem to have an effect on 
# any of the behaviors. Might need to inculde the stats for that in the results.
#note: although the function is 'mancova', just removing the covariate from the formula seems to give the correct results?
#note: there's good notation in the package for which test to use based on assumtions of normality and equal variance
# (e.g. pillai's trace v. wilk's lambda)
#2017

#best to use jmv package w/out covariate, results are identical to manova in base R, but all together in one neat table
mancova(data=behave.2017, deps=c('proportion.exposed','movements.min','bites.min',
                                 'total.dist.moved','courtship.min'),factors='Treatment')

#running this as manova in base R, not with jmv package
man.2017<-manova(cbind(proportion.exposed,movements.min,bites.min,total.dist.moved,courtship.min)~Treatment, data=behave.2017)
summary(man.2017)#no diff overall

# Look to see which differ, running univatiate contrasts
summary.aov(man.2017)#no diff in any of the univariate contrasts


#2018
mancova(data=behave.2018.t4.5, deps=c('proportion.exposed','movements.min','bites.min',
                                      'total.dist.moved','courtship.min'),factors='Treatment')

