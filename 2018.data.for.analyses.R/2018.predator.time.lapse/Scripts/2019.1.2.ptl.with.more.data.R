# Description: finished new data collection for ptl, wrangling and reanalyzing
# - going to wrangle raw(ish) data, in R script
# Author: George C Jarvis
# Date: Thu Jan 02 21:22:03 2020
# Notes: I went through and condensed the raw data down as much as I could:
#  A) I took out species, and condensed all observations down to a single row per reef
#  B) For analyses, data need to be in wide format, with separate categories for each photo class
#  C) For plotting, data need to be in short format (need column for "predator category")
# --------------

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
library(vegan)
library(multcomp)

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#raw data
pred<-read.csv(file = "Data/2019.1.2.ptl.raw.data.csv")
View(pred)
#removing all but first row of each observation per reef  (i.e. NA's)
pred.rm.na<-na.omit(pred)
head(pred.rm.na)

#exporting raw data, with each photo within a reef as an observation (n=10 photos minimum)

write.csv(pred.rm.na,"Data\\2019.1.2.predator.raw.data.csv", row.names = FALSE)
#write.csv(Your DataFrame,"Path where you'd like to export the DataFrame\\File Name.csv", row.names = FALSE)

View(pred.rm.na)#level of replication is photo per time-lapse per reef
#makes sense to me, because looking at proportion of photos from each time lapse
# that contains predators in different categories

#now want to average +/-, score of 4 (sublethal), and score of 5 (lethal),
# by Trial, Reef, Treatment.combo, and Treatment.t6 (for t6 comparison) 
#(dropping frame #, but will use to calc. SEM)

p<-with(pred.rm.na, aggregate(list(contained.pred.at.all.regardless.spp.,contianed.pred.score.4,contained.pred.score.5), 
                              list(Trial=Trial,Reef=Reef,Treatment.combo=Treatment.combo,Treatment.t6=Treatment.t6), mean))
write.csv(p,"Data\\2019.11.20.predatr.photos.averaged.csv", row.names = FALSE)


p
View(p) #I think it worked, but have to rename columns in new df

library(tidyverse)
p<-as_tibble(p)
p

#viewing column names to specify which ones I want to change
#columns correspond to following categories:

colnames(p)

#[presence] [4]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..1L..0L..0L.."  
#[sublethal] [5]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L.."  
#[lethal] [6]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L...1"

#renaming columns with base R
# Rename column names
names(p)[names(p) == "c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..1L..0L..0L.."] <- "Present"
names(p)[names(p) == "c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L.."] <- "Sublethal.Threat"
names(p)[names(p) == "c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L...1"] <- "Lethal.Threat"
p
View(p)

#data are now in wide format, 
# need to reshape to long format (p-long = "pl") for plotting

library(tidyr)

pl<-p %>% gather(Predator.class, Score, Present:Lethal.Threat)
View(pl)
#exporting in long format
write.csv(pl,"Data\\2019.11.20.predator.photos.long.format.csv", row.names = FALSE)


#new plot, no SEM though, not sure what SEM sciplot will use
pl$Treatment.combo<-ordered(pl$Treatment.combo, c("Low, Medium, High"))

bargraph.CI(x.factor = Treatment.combo, response = Score, group = Predator.class, legend=TRUE, main="predator presence, prelim",x.leg = 10, data = pl)

#looks similar to when I did it before

#will analyze with df in wide format (p), comaparing treatments with the same 
# achievable scores (i.e. 0-3 for all, 4 for MR and HR)

# changing df for analyses#####
#trying one variable at a time (presence -> sublethal -> lethal (really only applicable for T6))

#predator presence (applicable for all treatments)
p.presence<-with(pred.rm.na, aggregate((contained.pred.at.all.regardless.spp.), 
                                       list(Trial=Trial,Reef=Reef,Treatment.combo=Treatment.combo,Treatment.t6=Treatment.t6), mean))

#SEM for pred presence
#here, each photo within a time-lapse is the unit of replication
p.presence$pse<-with(pred.rm.na, aggregate((contained.pred.at.all.regardless.spp.), 
                                           list(Trial=Trial,Reef=Reef,Treatment.combo=Treatment.combo,Treatment.t6=Treatment.t6), function(x) sd(x)/sqrt(length(x))))[,5]
View(p.presence)

#sublethal presence (applicable for all but low risk)
p.sublethal<-with(pred.rm.na, aggregate((contianed.pred.score.4), 
                                        list(Trial=Trial,Reef=Reef,Treatment.combo=Treatment.combo,Treatment.t6=Treatment.t6), mean))

#SEM for pred presence
#here, each photo within a time-lapse is the unit of replication
p.sublethal$sse<-with(pred.rm.na, aggregate((contianed.pred.score.4), 
                                            list(Trial=Trial,Reef=Reef,Treatment.combo=Treatment.combo,Treatment.t6=Treatment.t6), function(x) sd(x)/sqrt(length(x))))[,5]
View(p.sublethal)

#lethal presence (applicable for high-risk caged, and high-risk uncaged)
p.lethal<-with(pred.rm.na, aggregate((contained.pred.score.5), 
                                     list(Trial=Trial,Reef=Reef,Treatment.combo=Treatment.combo,Treatment.t6=Treatment.t6), mean))

#SEM for pred presence
#here, each photo within a time-lapse is the unit of replication
p.lethal$lse<-with(pred.rm.na, aggregate((contained.pred.score.5), 
                                         list(Trial=Trial,Reef=Reef,Treatment.combo=Treatment.combo,Treatment.t6=Treatment.t6), function(x) sd(x)/sqrt(length(x))))[,5]
View(p.lethal)