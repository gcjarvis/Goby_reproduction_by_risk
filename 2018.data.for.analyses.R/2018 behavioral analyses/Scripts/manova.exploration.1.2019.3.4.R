# Description: manova for behavioral analyses, first take
# Author: George C Jarvis
# Date: Mon Mar 04 21:59:41 2019
# --------------

#NOTES: 1) I'm using density as a covariate in these analyses, and I've decided to use the density
#-- values recorded at the time of observation, and not those recorded for normal counts each week
#-- I think the measurements are not sure accurate if I use densities outside of the obs. period,
#-- and counted fish in the same way, regardless of whether I was doing an observation or a 
#-- density count for the reef (i.e. same methodology for both)

#2) That being said, I did have to use density data for 2017 trials from my main density data
#-- because I didn't record densities prior to observations in 2017

#3) Need to go back to raw data sheet to reconcile some things
#--a) not sure if blank cells in behavioral datasheet for bites + distance 
#-- moved are 0's or missing observations
#--b) not sure if the distance moved column is based off of cm or mm measurements
#-- I think I made a note of it in the 2018 raw datasheets at one point
#-- and I feel like it was labeled in the 2017 Excel file for observation template
#-- confirmed: 2017 datasheet says "distance.moved.(mm)" BUT
#-- HB recorded some fish moving 1250mm, which = 4 slate lengths


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

#2019.3.12 data with 2017 densities:
behave<-read.csv("Data/behavior.data.2019.3.14.csv") #doesn't include unnecessary NA's
#behave<-read.csv("Data/behavior.data.2019.3.12.csv") #includes data on gravidity, and tags of chaser/chasee
#data read in for packet of analyses: behave<-read.csv("Data/behavior.data.3.4.19.added.zeros.for.blanks.bites.csv")
#ordering treatments for better graphical visualization
#behave$Treatment<-ordered(behave$Treatment,levels=c("Low","Medium","High","Control"))

#subsetting data, only 2017
behave.2017<-behave[behave$Trial<4, ]
#behave.2017<-behave[behave$Trial<4, c("Year","Trial", "Reef", "Week","Days.after.deployment","Observer","Density","Goby.ID", "Treatment", "Density","Gravidity","Nearest.TOL","Starting.location",
#                                         "Starting.position","Total.time.exposed","proportion.exposed","Total.time.hidden","proportion.hidden","num.movement.swims","num.bites","num.times.chased",
#                                        "chased.by.size.tag","num.times.chaser","chaser.of.size.tag","area.used.cm","Distance.moved.review.this.column","corrected.distance.moved")]
#2018 data, only trials 4 and 5
behave.2018.t4.5<-behave[c((behave$Trial>3) & (behave$Trial<6)), ]

#behave.2018.t4.5<-behave[c((behave$Trial>3) & (behave$Trial<6)), c("Year","Trial", "Reef", "Week","Days.after.deployment","Observer","Density","Goby.ID", "Treatment", "Density","Gravidity","Nearest.TOL","Starting.location",
#                                         "Starting.position","Total.time.exposed","proportion.exposed","Total.time.hidden","proportion.hidden","num.movement.swims","num.bites","num.times.chased",
#                                         "chased.by.size.tag","num.times.chaser","chaser.of.size.tag","area.used.cm","Distance.moved.review.this.column","corrected.distance.moved")]

#2018 data, only trial 6
behave.2018.t6<-behave[behave$Trial>5, ]
#behave.2018.t6<-behave[behave$Trial>5, c("Year","Trial", "Reef", "Week","Days.after.deployment","Observer","Density","Goby.ID", "Treatment", "Density","Gravidity","Nearest.TOL","Starting.location",
#                                                                   "Starting.position","Total.time.exposed","proportion.exposed","Total.time.hidden","proportion.hidden","num.movement.swims","num.bites","num.times.chased",
#                                                                   "chased.by.size.tag","num.times.chaser","chaser.of.size.tag","area.used.cm","Distance.moved.review.this.column","corrected.distance.moved")]
#making behave = every trial df
behave<-behave.2017
behave<-behave.2018.t4.5
behave<-behave.2018.t6

View(behave)

#complete cases only
behave<-behave[complete.cases(behave),]


#differences in behavior among treatments?
mod1<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                   corrected.distance.moved)~Treatment,data=behave)
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
summary(mod1)

#see which differ among treatments
summary.aov(mod1)
#looks like differences are driven by exposure time and distance moved,
#although I am less confident in the distance moved metric (cm vs. mm
# discrepency in dataset)

#pairwaise comparisons to see which treatments drove trends
#1) All data combined
#low vs. medium = NS
mod1a<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                   Distance.moved.review.this.column)~Treatment,data=behave,
            subset = Treatment %in% c("Low","Medium"))
summary(mod1a)

#low vs. high = sig.
mod1b<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                    Distance.moved.review.this.column)~Treatment,data=behave,
              subset = Treatment %in% c("Low","High"))
summary(mod1b)

#high vs. medium = sig.
mod1c<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                    Distance.moved.review.this.column)~Treatment,data=behave,
              subset = Treatment %in% c("Medium","High"))
summary(mod1c)

#high vs. control = NS
mod1d<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                    Distance.moved.review.this.column)~Treatment,data=behave,
              subset = Treatment %in% c("High","Control"))
summary(mod1d)

#...can do more comparisons, but next move should be to break it 
#--down by trial to see if trends hold, or if we run into ss issues (trial 6)


#differences in behavior among trials?
mod2<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                 Distance.moved.review.this.column)~as.factor(Trial),data=behave)
summary(mod2)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
boxplot(Total.time.exposed*num.movement.swims*num.bites*
        Distance.moved.review.this.column~Trial,data=behave)

#all measured behaviors differed among trials
summary.aov(mod2)
#pairwise comparisons to see which trials differed (2017 first)
#1v2- not diff
mod1a<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                    Distance.moved.review.this.column)~Treatment,data=behave,
              subset = Trial %in% c("1","2"))
summary(mod1a)

#1v3-not diff
mod1b<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                    Distance.moved.review.this.column)~Treatment,data=behave,
              subset = Trial %in% c("1","3"))
summary(mod1b)

#2v3-not diff
mod1c<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                    Distance.moved.review.this.column)~Treatment,data=behave,
              subset = Trial %in% c("2","3"))
summary(mod1c)

#high vs. control = NS
#mod1d<-manova(cbind(Total.time.exposed,num.movement.swims,num.bites,
                    #Distance.moved.review.this.column)~Treatment,data=behave,
              #subset = Trial %in% c("High","Control"))
#summary(mod1d)

##MANCOVA from Youtube video####

#testing assumption of parallel slopes is true 
#if interaction term is insignificant, then assumptions of parallel slopes is true
outcome<-cbind(behave$Total.time.exposed,behave$num.movement.swims,behave$num.bites,
              behave$Distance.moved.review.this.column)
outcome
model<-manova(outcome~Treatment+Density+Treatment*Density, data=behave)
summary(model,test = "Wilks", type="III")
summary(model, test="Pillai", type="III")
summary(model,test="Hotelling", type="III")
summary(model, type="Roy", type="III")
#for all trials combined, assumptions of parallel slopes is true

#looking at means of outcomes in each group
#uses psych package
#crazy output for entire combined dataset, but might be nice when just looking
#at condensed datasets among trials
describeBy(behave,behave$Treatment)

#now want to see if adjusted means for each behavior differ among treatments

#calsulate adjusted means with aov() function and effect() function
#from "effect" package
#have to go through and calculte adjusted means for each factor
factor(behave$Treatment)
#calculating adjusted means
#significant results show thatfactor has a sig. effect on outcome variable
#ex. mod.expose shows that density and treatment sig. affected exposure time
mod.expose<-aov(Total.time.exposed~Density+Treatment,data=behave)
summary(mod.expose, type="III")

adjmean.expose<-effect("Treatment",mod.expose,se=TRUE,xlevels=4)#I thinmk I should keep levels at 4
#although in her video, woman keeps levels at 2
summary(adjmean.expose)
adjmean.expose$se

#adjusted means for exposure time, with levels as 2
#(control, high, low, medium)
#(238.6678, 214.2169, 249.6051, 233.4091)

#same but with levels as 4
#(165.8895 196.2357 233.0024 217.0068)

#looking at next factor of movement swims

#no diff in # of movement swims by treatment, but density had an effect (all data)
mod.movements<-aov(num.movement.swims~Density+Treatment,data=behave)
summary(mod.movements, type="III")

adjmean.movements<-effect("Treatment",mod.movements,se=TRUE,xlevels=4)#I think I should keep levels at 4
#although in her video, woman keeps levels at 2
summary(adjmean.movements)
adjmean.movements$se

#adj. means for movement swims
#(5.431169 6.032638 6.788719 5.612326)

#bites
mod.bites<-aov(num.bites~Density+Treatment,data=behave)
summary(mod.bites, type="III")
#all data: bites differ by density, not by treatment

adjmean.bites<-effect("Treatment",mod.bites,se=TRUE,xlevels=4)#I think I should keep levels at 4
#although in her video, woman keeps levels at 2
summary(adjmean.bites)
adjmean.bites$se

#adj. means for bites
#(2.298245 4.033894 3.524353 2.887672)

#distance moved
mod.dist<-aov(num.bites~Distance.moved.review.this.column+Treatment,data=behave)
summary(mod.dist, type="III")
#all data: no effect of treatment or density on distance moved

adjmean.dist<-effect("Treatment",mod.dist,se=TRUE,xlevels=4)#I think I should keep levels at 4
#although in her video, woman keeps levels at 2
summary(adjmean.dist)
summary(adjmean.dist$se)

#adj.means dist.moved
#(1.900677 3.993359 3.306725 3.322092)

#MANCOVA with JMV package######
#seems like this could be the way to go, even though it might not satisfy all assumptions
#library(jmvcore)
library(jmv)
#mancova(data=behave, deps=c('Total.time.exposed','num.movement.swims','num.bites',
#                            'Distance.moved.review.this.column'),factors='Treatment',covs = 'Density')

#2017
#NOTE: I didn't record densities during my 2017 observations
#added in density data in this manner:
#--1) if I recorded density for the day of the behavioral obs. then I used that density value
#--2) if I did not record density for the obs. day, then I averaged the densities taken
#--   from the two days around the observation day
#--3) I recorded the method from where I got density in the .csv file as "Method.for.den
#--   I only had to calculate average density for 5 of the reefs in 2017 (all in trial 2)
#--from observations only, which would allow me to enter these data into 2017 densities

#went back and adjusted the distance moved metric becuase I was uncertain 
#--as to whether the units were off. Called the new variable "corrected.distance.moved"

#also added in number of times chased and number of times chaser

#original MANCOVA
mancova(data=behave, deps=c('Total.time.exposed','num.movement.swims','num.bites',
                            'Distance.moved.review.this.column','num.times.chased',
                            'num.times.chaser'),factors='Treatment',covs='Density')

#MANCOVA with corrected distances WITHOUT original distance moved
#went with this analysis for the results that I reported to Mark
mancova(data=behave1, deps=c('Total.time.exposed','num.movement.swims','num.bites',
                            'corrected.distance.moved','num.times.chased',
                            'num.times.chaser'),factors='Treatment',covs='Density')

#look at ss by treatment:

#READ THIS RE: SS! for my thesis, I counted the number of individuals observed
#--for each treatment within each experiment as the sample size
#--Sample sizes were really big for exp 1 and 2, but small for exp 3

library(plyr)

#count(behave, "Treatment")
#2017 w/ complete cases:
#n = (high-54, low-55, medium-53)

#how to get sample size by trial and treatment with xtabs:
#xtabs( ~ Treatment + Reef + Trial, behave) 
#NOTE: I observed at least 3 individuals per reef per treatment
#---functionally, I think MANCOVA is averaging behaviors for each 
#---treatment, which is okay, despite low ss in trial 6?

#2017: 6 observations per treatment per trial, n=18 for each
#2018 t4 and 5: t4: l: 6, m: 5, h: 6, t5: l: 6, m: 6, h: 6.
#-----totals: l:12 , m:11 , h:12
#2018 t6: l:1 , m:2 , h:1 , c:2 

#MANCOVA with corrected distances WITH original distance moved, doesn't seem to work
#mancova(data=behave, deps=c('Total.time.exposed','num.movement.swims','num.bites',
#                            'Distance.moved.review.this.column','corrected.distance.moved','num.times.chased',
#                           'num.times.chaser'),factors='Treatment',covs='Density')

#2018.t.4.5
mancova(data=behave, deps=c('Total.time.exposed','num.movement.swims','num.bites',
                            'Distance.moved.review.this.column'),factors='Treatment',covs = 'Density')
#did gobies behave as expected based on treatment?
#looking at behaviors that differed within MANCOVA: time exposed, and distance moved
# exposure time
tapply(behave$Total.time.exposed, behave$Treatment, mean)
summary(mod1<-aov(Total.time.exposed~Treatment*Density,data=behave))
TukeyHSD(mod1,"Treatment",ordered = TRUE)
plot(TukeyHSD(mod1,"Treatment"))

mod1<-aov(Total.time.exposed~Treatment*Density,data=behave)
plot(mod1)

#n per treatment with plyr package
count(behave, "Treatment")
#n = (high-76, low-98, medium-102)

# distance moved
#-- need to deal with NA values first
behave.a<-behave[!is.na(behave$Distance.moved.review.this.column),]
tapply(behave.a$Distance.moved.review.this.column, behave.a$Treatment, mean)
summary(mod1<-aov(Distance.moved.review.this.column~Treatment*Density,data=behave.a))
TukeyHSD(mod1,"Treatment",ordered = TRUE)
plot(TukeyHSD(mod1,"Treatment"))

mod1<-aov(Distance.moved.review.this.column~Treatment*Density,data=behave.a)
plot(mod1)

#2018.t6
mancova(data=behave, deps=c('Total.time.exposed','num.movement.swims','num.bites',
                            'Distance.moved.review.this.column'),factors='Treatment',covs = 'Density')

#n per treatment with plyr package
count(behave, "Treatment")
#n = (control-6,high-3, low-3, medium-6)