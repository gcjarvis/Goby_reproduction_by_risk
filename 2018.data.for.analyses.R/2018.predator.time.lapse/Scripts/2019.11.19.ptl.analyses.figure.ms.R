# Description: PTL analyses and figure for MS
# Author: George C Jarvis
# Date: Tue Nov 19 09:23:42 2019
# Notes: Can only compare HR-caged vs. HR-uncaged within trial 6, separately
#       Will have to do that after
#       objectives:
#       1) wrangle raw data into better df's 
#       2) figure out the stats for PERMANOVA, and where it might have gone wrong
#       3) plotting, and how to make grayscale
#
#   NOTE: I think I might have to break the data down into 4 separate analyses:
#       1) presence/absence among all treatments
#       2) sublethal predation between med and HR (combined caged and uncaged)
#       3) lethal predation is going to be present in HR, so different than 0 = statistically sig?
#       4) high-risk caged vs. high-risk uncaged?
# 
#
# Analyses: I think wide-format is best way to have df, plotting will be long format
#       - ****need to figure out how to get SEM, may need to break down each predator classification one at a time
#         - instead of doing it all at once with aggregate function
# --------------

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

#raw data
pred<-read.csv(file = "Data/2019.11.19.ptl.updated.data.csv")
#removing Count.zero and Score.zero (all NA's)
pred.rm<-subset(pred, select = -c(Count, Score))
pred.rm
#removing all but first row of each observation per reef  (i.e. NA's)
# (i.e. where I calculated counted +/- preds for different positional bins)
pred.rm.na<-na.omit(pred.rm)
head(pred.rm.na)

#exporting raw data, with each photo as an observation

write.csv(pred.rm.na,"Data\\predator.photos.raw.csv", row.names = FALSE)
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

#combining df's, not sure how I will plot this, may have to make separat plots and overlay them
#THIS DOESN'T SEEM TO BE WORKING, DATA LOOK OFF TO ME
p.combined<-left_join(p.presence, p.sublethal,p.lethal, by = c("Trial","Reef","Treatment.combo","Treatment.t6"))
View(p.combined)
#lost lethal x and se
#trying full join instead
p.full.join<-full_join(p.presence, p.sublethal,p.lethal, by = c("Trial","Reef","Treatment.combo","Treatment.t6"))

View(p.full.join)

left_join(a, b, by = "x1")

#selecting trial 6 data only



