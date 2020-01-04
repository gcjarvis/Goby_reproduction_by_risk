# Description: finished new data collection for ptl, wrangling and reanalyzing
# - going to wrangle raw(ish) data, in R script
# Author: George C Jarvis
# Date: Thu Jan 02 21:22:03 2020
# Notes: I went through and condensed the raw data down as much as I could:
#  A) I took out species, and condensed all observations down to a single row per reef
#  B) For analyses, data need to be in wide format, with separate categories for each photo class
#  C) For plotting, data need to be in short format (need column for "predator category")
# Lastly, it's important to note how these are being coded for figures:
# if a single predator was seen in a photo, it was marked as a 1 for containing predator
# That's regardless of where the predator was seen in the frame
# I mention this because if a predator were seen in a photo and at the edge of a sm. cage
# - then it would count for both present, and sublethal risk
# I did not, however, count lethal predators as being both sublethal and lethal (i.e. they
# - were only counted as one or the other)
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

#type III ANOVA
options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#importing data####
#raw data
pred<-read.csv(file = "Data/2019.1.2.ptl.raw.data.csv")
View(pred)
#removing NA's
pred.rm.na<-na.omit(pred)
head(pred.rm.na)

#adding two columns for year, one for year (int), and one for year.fact (fact) to df;
# - as a proxy for tagging procedure, where trials 1-3 = 2017, and 4-6 = 2018
pred.rm.na$Year <- ifelse(pred.rm.na$Trial <=3, 2017, 2018)
#want it as a factor? Going to make a variable with it as a factor, run the model, and see if I get different results
pred.rm.na$Year.fact<- as.factor(pred.rm.na$Year)
View(pred.rm.na)

#exporting raw data, with each photo within a reef as an observation (n=10 photos minimum)

# These are the data in wide format
write.csv(pred.rm.na,"Data\\2019.1.2.predator.raw.data.csv", row.names = FALSE)
#write.csv(Your DataFrame,"Path where you'd like to export the DataFrame\\File Name.csv", row.names = FALSE)

View(pred.rm.na)#level of observation is photo per time-lapse per reef, but replicate will be avg. proportion per reef

#importing egg count data for avg. number of gobies inhabiting each reef
# Will add this to the data when I wrangle them into wide format
repro<-read.csv(file= "Data/egg.counts.2019.12.23.csv")
#making Year.fact a factor
repro$Year.fact<-as.factor(repro$Year.fact)

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
                         T6.comparison=T6.comparison,Year=Year,Year.fact=Year.fact), mean))
View(p)

#Column names are not correct, need to change them

#viewing column names to specify which ones I want to change
#columns correspond to following categories:

colnames(p)

#[presence] [7]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..1L..0L..0L.."  
#[sublethal] [8]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L.."  
#[lethal] [9]"c.0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L..0L...1"

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

pr<-left_join(p,repro,by=c("Trial","Reef","Treatment","T6.comparison","Year","Year.fact"))
View(pr)

#want to drop unnecessary columns from repro df, i.e. everything but avg.inhab
# they are columns 8:12 and 14:16
pr1<-pr[-c(10:14,16)]
View(pr1)

# this looks good to me, checked with repro data to make sure avg.inhab values were correct

#exporting data in wide format
write.csv(pr1,"Data\\2019.1.3.predator.data.wide.format.with.avg.inhab.csv", row.names = FALSE)

#long format for plotting#### 
# p-long = "pl"

#for all trials
pr1<-as_tibble(pr1)
pr1

pl<-pr1 %>% gather(Predator.class, Score, Present:Lethal.Threat)
View(pl)
#exporting in long format
write.csv(pl,"Data\\2019.1.3.predator.data.long.format.with.avg.inhab.csv", row.names = FALSE)

#separate df named "HR.comp" to compare responses between HR-caged and -uncaged treatments
# Will not include year or trial in the model for analyses
# I only did time lapses on high-risk caged and uncaged treatments for T6

HR.comp<-pr1[pr1$Trial==6,]
View(HR.comp)
#exporting data
write.csv(HR.comp,"Data\\2019.1.2.predator.raw.data.high.risk.comparison.csv", row.names = FALSE)

#NOTE: might have to go back and just do Trials < 6 to compare all others
# that depends on whether there is a trial effect

#doing exactly this now:
#data frame for trials 1-5 in wide format
t1.5.w<-pr1[pr1$Trial<6,]
View(t1.5.w)
#exporting data
write.csv(t1.5.w,"Data\\2019.1.2.predator.raw.data.trials.1.5.wide.csv", row.names = FALSE)

# ^ this df will work for analyss for presence/absence, but not for sublethal or lethal
#need to make df's to subset by treatment

#going to drop low from the model and see if there are still differences
ptl.sub<-subset(t1.5.w,Treatment!="Low")
View(ptl.sub)
levels(ptl.sub$Treatment)

#no formal analysis will be done for lethal, just make sure to plot it and
# - include in the figure legend that no statistical test was run for high-risk caged
# - in trials 1-5

#ask M. Steele: if there are differences between caged and uncaged trts
# - (e.g. presence of preds) do I have to remove that trial from the 
# - overall analyses for all trials, and run trials 1-5 and trial 6 as 
# - two separate anlayses? I think yes (see df "t1.5.l")

# will see if there are statistical differences among trials

#for trials 1-5 only (l, m, h (caged) only)
trial.1.5.long<-pl[pl$Trial<6,]
View(trial.1.5.long)

#for trial 6 only (HR caged and uncaged)
HR.long<-pl[pl$Trial==6,]
View(HR.long)

#calculating sample sizes and number of photos that were used to calc. proportions per time lapse####
# average number of photos that went into each proportion per reef per treatment

# a) including uncaged
#had to do it in excel pivot table...

#avg.number.of.photos
#high	11.48 (12)
#med	11.6 (12) 
#low	12.1 (12)
#uncaged	13.5 (14)

# b) with HR caged and uncaged combined

#avg.number.of.photos
#high	11.8 (12)
#med	11.6 (12) 
#low	12.1 (12)

#arithmetic means using wide format df (pr1), mainly for ss of treatments
library(FSA)

#with uncaged, only compared caged and uncaged high-risk in final trial
Summarize(Present~T6.comparison,
          data=HR.comp,
          digits=3)
# n:
# H - 4
# U - 4

#with uncaged and caged combined, not thinking this is correct, see below
Summarize(Present~Treatment,
          data=pr1,
          digits=3)
# n:
# L - 21
# M - 20
# H - 29

#with HR caged only, using t1.4.w df, not including any data from t6
Summarize(Present~Treatment,
          data=t1.5.w,
          digits=3)
# n:
# L - 21
# M - 20
# H - 21

#using df "pr1" --> predator data in wide format


#all trials
#going to evaluate predator activity with nested factor of trial - bring in avg.inhab? 
# It might be compelling to bring that in, considering there preds might be attracted to reefs 
# - with more fish on them?

#comparison between high-risk caged and uncaged (trial 6 only)

#analyses####
#NOTE: will need to test results among trials to see if there are any sig. differences based
# - on trial alone, which might not allow me to group data from all trials together

# 1) will compare high-risk caged and uncaged treatments

# a) proportion of photos that contained predators at all (comparable among all trt's)
#full model, including covariate of avg.inhab
# using reduced df for comparison --> "HR.comp"

modp<-lm(Present~T6.comparison*avg.inhab, data=HR.comp)
qqnorm(resid(modp))
qqline(resid(modp))
summary(modp)
anova(modp) # interactive effects of avg.inhab and cage on proportion of photos
# - where predators were present

#my instinct is to not combine HR caged and uncaged for this analysis, but maybe for
# - others where there is no effect of cages

#see plotting for an idea of what was going on with interaction
#more goboes = higher prop. of photos with preds, but only when HR reefs were caged

#seems to be an effect of avg. inhab on proportion of photos with predators
#plotting relationship between avg. inhab and proportion of photos with predators
plot(HR.comp$avg.inhab, HR.comp$Present)
abline(lm(HR.comp$Present~HR.comp$avg.inhab*HR.comp$T6.comparison))
# I think barplot is the best comparison for this, 
# - but will check lineplot too (see plotting section)

#b) proportion of photos with predators as sublethal threats
mods<-lm(Sublethal.Threat~T6.comparison*avg.inhab, data=HR.comp)
qqnorm(resid(mods))
qqline(resid(mods))
summary(mods)
anova(mods)

#bi) removing interactions with covariate from the model
mods.1<-lm(Sublethal.Threat~T6.comparison+avg.inhab, data=HR.comp)
qqnorm(resid(mods.1))
qqline(resid(mods.1))
summary(mods.1)
anova(mods.1)

#bii) removing covariate from the model completely
mods.2<-lm(Sublethal.Threat~T6.comparison, data=HR.comp)
qqnorm(resid(mods.2))
qqline(resid(mods.2))
summary(mods.2)
anova(mods.2)

#c) proportion of photos containing lethal threat
modl<-lm(Lethal.Threat~T6.comparison*avg.inhab, data=HR.comp)
qqnorm(resid(modl))
qqline(resid(modl))
summary(modl)
anova(modl)
#strange that SS, MS, anf F-value are 0; checking df
view(HR.comp)

#ci) reducing model by removing interaction
modl.1<-lm(Lethal.Threat~T6.comparison+avg.inhab, data=HR.comp)
qqnorm(resid(modl.1))
qqline(resid(modl.1))
summary(modl.1)
anova(modl.1)

#cii) removing covariate completely
modl.2<-lm(Lethal.Threat~T6.comparison, data=HR.comp)
qqnorm(resid(modl.2))
qqline(resid(modl.2))
summary(modl.2)
anova(modl.2)
#want to plot this out and see what's up with model

# 1) proportion of photos that contained predators at all (comparable among all trt's)
#NOTE: have to make a new df that is only trials 1-5, because prop.present varied
#between caged and uncaged treatments in T6

#I'm going to keep the models that I ran initially, BUT going to run
# - "better" model with trial 1-5 only, where all treatments were present

#analyses run with new df's for trials 1-5 only (df = "t1.5.w")
#BIG NOTE: these are the models that I put in the MS draft
# - that I was working on on 2019.1.3-

#ai) full model
modpi<-lme(Present~Treatment*Year.fact*avg.inhab,
           random=~1|Trial,t1.5.w,method="REML")
hist(resid(modpi))
qqnorm(resid(modpi))
qqline(resid(modpi))
summary(modpi)
anova(modpi, type='marginal')

#aii) reduced, no 3-way
modpii<-lme(Present~Treatment*Year.fact+(Treatment*avg.inhab)+(Year.fact*avg.inhab)+
          avg.inhab,random=~1|Trial,t1.5.w,method="REML")
hist(resid(modpii))
qqnorm(resid(modpii))
qqline(resid(modpii))
summary(modpii)
anova(modpii, type='marginal')

#aiii) reduced no 2-ways with covariate
modpiii<-lme(Present~Treatment*Year.fact+
              avg.inhab,random=~1|Trial,t1.5.w,method="REML")
hist(resid(modpiii))
qqnorm(resid(modpiii))
qqline(resid(modpiii))
summary(modpiii)
anova(modpiii, type='marginal')

#aiv) reducing further, no covariate at all
modpiv<-lme(Present~Treatment*Year.fact,random=~1|Trial,t1.5.w,method="REML")
hist(resid(modpiv))
qqnorm(resid(modpiv))
qqline(resid(modpiv))
summary(modpiv)
anova(modpiv, type='marginal')

#models run INITIALLY, but I think they are wrong 
# - b/c they assume that hr caged and uncaged are the same
#a) full model, including covariate of avg.inhab
modpa<-lme(Present~Treatment*Year.fact*avg.inhab,random=~1|Trial,pr1,method="REML")
qqnorm(resid(modpa))
qqline(resid(modpa))
summary(modpa)
anova(modpa, type='marginal')

#ai) removing three-way interaction with covariate
modpa.1<-lme(Present~Treatment*Year.fact+(Treatment*avg.inhab)+
            (Year.fact*avg.inhab)+avg.inhab,random=~1|Trial,
            pr1,method="REML")
qqnorm(resid(modpa.1))
qqline(resid(modpa.1))
summary(modpa.1)
anova(modpa.1, type='marginal')

#aii) removing interactions with covariate
modpa.2<-lme(Present~Treatment*Year.fact+avg.inhab,random=~1|Trial,
             pr1,method="REML")
qqnorm(resid(modpa.2))
qqline(resid(modpa.2))
summary(modpa.2)
anova(modpa.2, type='marginal')

#aiii) removing covariate altogether
modpa.3<-lme(Present~Treatment*Year.fact,random=~1|Trial,
             pr1,method="REML")
qqnorm(resid(modpa.3))
qqline(resid(modpa.3))
summary(modpa.3)
anova(modpa.3, type='marginal')
ranef(modpa.3)

#testing for differences among trials, including as a fixed effect
modpa.3t<-aov(Present~Treatment+Year.fact+avg.inhab+Trial, data=pr1)
hist(resid(modpa.3t))
qqnorm(resid(modpa.3t))
qqline(resid(modpa.3t))
summary(modpa.3t)
anova(modpa.3t)

# I tested for all interactions as fixed effects, and it seems that there
# - weren't any differences in predator presence among trials
# I'd say run it with the same model as egg counts

# 2) proportion of photos that contained sublethal predators
# (comparable among medium- and high-risk treatments only)
#not quite sure if the subset df that I made only contians comparable treatments
# - i.e., it liiks like low is still included, but will try and model it for a test

#NOTE: these are the data that I put in the table for the MS draft (2019.1.3)
#a) full model, including covariate of avg.inhab
modsi<-lme(Sublethal.Threat~Treatment*Year.fact*avg.inhab,
           random=~1|Trial,ptl.sub,method="REML")
hist(resid(modsi))
qqnorm(resid(modsi))
qqline(resid(modsi))
summary(modsi)
anova(modsi, type='marginal')
fixef(modsi)

#checking this out further, seems to be working correctly
Summarize(Sublethal.Threat~Treatment,
          data=ptl.sub,
          digits=3)

#ai) removing three-way interaction with covariate
modsii<-lme(Sublethal.Threat~Treatment*Year.fact+(Treatment*avg.inhab)+
            (Year.fact*avg.inhab)+avg.inhab,random=~1|Trial,
            ptl.sub,method="REML")
hist(resid(modsii))
qqnorm(resid(modsii))
qqline(resid(modsii))
summary(modsii)
anova(modsii, type='marginal')

#aii) removing interactions with covariate
modsiii<-lme(Sublethal.Threat~Treatment*Year.fact+avg.inhab,random=~1|Trial,
            ptl.sub,method="REML")
hist(resid(modsiii))
qqnorm(resid(modsiii))
qqline(resid(modsiii))
summary(modsiii)
anova(modsiii, type='marginal')

#aiii) removing covariate altogether
modsiv<-lme(Sublethal.Threat~Treatment*Year.fact,random=~1|Trial,
             ptl.sub,method="REML")
hist(resid(modsiv))
qqnorm(resid(modsiv))
qqline(resid(modsiv))
summary(modsiv)
anova(modsiv, type='marginal')

#testing for differences among trials, including as a fixed effect

#interactions
modsv<-aov(Present~Treatment*Year.fact*avg.inhab*Trial, data=ptl.sub)
hist(resid(modsv))
qqnorm(resid(modsv))
qqline(resid(modsv))
summary(modsv)
anova(modsv) #nothing significant

#no interactions
modsvi<-aov(Present~Treatment+Year.fact+avg.inhab+Trial, data=ptl.sub)
hist(resid(modsvi))
qqnorm(resid(modsvi))
qqline(resid(modsvi))
summary(modsvi)
anova(modsvi)

#doesn't seem to be an effect of trial on the prop of photos that contained sublethal predators


#plotting, using data in long format (trials 1-5: pl, trial6 6: HR.comp)####

#A) diagnostic plots to see how caged and uncaged treatments compared in trial 6

#1) predator presence
bargraph.CI(x.factor = avg.inhab, response = Present, group = T6.comparison, 
            legend=TRUE, main="predator presence, t6 comp", data = HR.comp)

lineplot.CI(avg.inhab,Present,group=T6.comparison,legend = TRUE,
            main="predator presence, t6 comp", xlab="avg.inhab",
            ylab="proportion of photos predator present", data=HR.comp)



bargraph.CI(x.factor = Treatment.combo, response = Score, group = Predator.class, legend=TRUE, main="predator presence, prelim",x.leg = 10, data = pl)

#2) sublethal threat
#not really needed, no differences statistically, but going to combine lethal and sublethal
#into one frame

#first, subset df to only include lethal and sublethal from HR.long
HR.long.sub<-subset(HR.long,Predator.class!="Present")
HR.long.sub<-subset(HR.long.sub,T6.comparison!=c("Low","Medium"))
View(HR.long.sub)


#plotting
bargraph.CI(x.factor = T6.comparison, response = Score, group = Predator.class,
            legend=TRUE, main="predator lethal + sublethal, t6 comp", data = HR.long.sub)


#3) lethal threat
#no differences, but want to see what's going on with weird F-value

bargraph.CI(x.factor = T6.comparison, response = Lethal.Threat,
            legend=TRUE, main="predator lethal, t6 comp", data = HR.comp)


#B) diagnostic plots for trials 1-5 (all treatments)
#will have to subset data to make comparisons
# just do boxplots to visualize? Can't do a formal comparison for lethal threat
# - because it was only possible in HR-caged
pl$Treatment.combo<-ordered(pl$Treatment.combo, c("Low, Medium, High"))

#looks similar to when I did it before

#will analyze with df in wide format (p), comaparing treatments with the same 
# achievable scores (i.e. 0-3 for all, 4 for MR and HR)

#plotting figures, will likely be my final plots in ggplot
#doing with barplot now

#trials 1-5, df = "trial.1.5.long"
trial.1.5.long$Treatment.ordered<-ordered(trial.1.5.long$Treatment,c("Low","Medium","High"))
trial.1.5.long$Predator.class.ordered<-ordered(trial.1.5.long$Predator.class,c("Present","Sublethal.Threat","Lethal.Threat"))
bargraph.CI(x.factor = Treatment.ordered, response = Score, group = Predator.class.ordered,
            legend=TRUE, xlab="Treatment", ylab= "proportion of photos", main="predator activity, trials 1-5", data = trial.1.5.long)

View(trial.1.5.long)

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