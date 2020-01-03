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

#packages####
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
library(tidyverse)
library(tidyr)
library(qpcR)#for cbind.na

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#importing data####
#raw data
pred<-read.csv(file = "Data/2019.1.2.ptl.raw.data.csv")
View(pred)
#removing all but first row of each observation per reef  (i.e. NA's)
pred.rm.na<-na.omit(pred)
head(pred.rm.na)

#exporting raw data, with each photo within a reef as an observation (n=10 photos minimum)

# These are the data in wide format
write.csv(pred.rm.na,"Data\\2019.1.2.predator.raw.data.csv", row.names = FALSE)
#write.csv(Your DataFrame,"Path where you'd like to export the DataFrame\\File Name.csv", row.names = FALSE)

View(pred.rm.na)#level of replication is photo per time-lapse per reef
#makes sense to me, because looking at proportion of photos from each time lapse
# that contains predators in different categories

#importing egg count data for avg. number of gobies inhabiting each reef
# Will add this to the data when I wrangle them into wide format
repro<-read.csv(file= "Data/egg.counts.2019.12.23.csv")

#data wrangling####

#wide format (for analyses)####

# want the proportion of photos per reef +/- predators, 
# - predators close enough to be a sublethal threat, 
# - and predators close enough to be a lethal threat

# Grouping by Trial, Reef, Treatment.combo, and Treatment.t6 (for t6 comparison) 
# Read number is the level of observation
# Each reef is a replicate, so that is what will be used to calculate means and SEM

p<-with(pred.rm.na, aggregate(list(contained.pred,sublethal.threat,lethal.threat), 
                    list(Trial=Trial,Reef=Reef,Treatment=Treatment,
                         T6.comparison=T6.comparison), mean))
View(p)

#Column names are not correct, need to change them

#viewing column names to specify which ones I want to change
#columns correspond to following categories:

colnames(p)

#[presence] [5]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..1L..0L..0L.."  
#[sublethal] [6]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L.."  
#[lethal] [7]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L...1"

#renaming columns with base R
# Rename column names
names(p)[names(p) == "c.0L..0L..1L..0L..1L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L.."] <- "Present"
names(p)[names(p) == "c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L.."] <- "Sublethal.Threat"
names(p)[names(p) == "c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L...1"] <- "Lethal.Threat"
p
View(p)

#want to bring in avg. numer of gobies inhabiting each reef as well, adding it to "p" df
#bringing in egg count data, so I can assign avg.inhab variable to the corresponding rows in my p df
# - There are uneven #'s of rows between the two df's, so I'll see what happens with that

pr<-left_join(p,repro,by=c("Trial","Reef","Treatment","T6.comparison"))
View(pr)

#want to drop unnecessary columns from repro df, i.e. everything but avg.inhab
# they are columns 8:12 and 14:16
pr1<-pr[-c(8:12,14:16)]
View(pr1)

# this looks good to me, checked with repro data to make sure avg.inhab values were correct

#exporting data in wide format
write.csv(pr1,"Data\\2019.1.3.predator.data.wide.format.with.avg.inhab.csv", row.names = FALSE)

# pr1 contains the data to be used in the analyses for ptl

#long format for plotting#### 
# p-long = "pl"

pr1<-as_tibble(pr1)
pr1

pl<-pr1 %>% gather(Predator.class, Score, Present:Lethal.Threat)
View(pl)
#exporting in long format
write.csv(pl,"Data\\2019.1.3.predator.data.long.format.with avg.inhab.csv", row.names = FALSE)

#analyses####

#using df "pr1" --> predator data in wide format


#all trials
#going to evaluate predator activity with nested factor of trial - bring in avg.inhab? 
# It might be compelling to bring that in, considering there preds might be attracted to reefs 
# - with more fish on them?

#comparison between high-risk caged and uncaged (trial 6 only)

#plotting, using pl df, in long format####

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