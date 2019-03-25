# Description: predator time lapse, initial analyses  
# Author: George C Jarvis
# Date: Sat Mar 02 14:52:22 2019
# --------------

#objectives:
# 1. determine effectiveness of caging treatments based on average position scores of predators
# 2. determine makeup of predator assemblages
# 3. determine if there is any difference in the number of fish seen voer time
#-- this will help me figure out if preds. were afraid of my presence immediately after
#-- deploying the camera

#general notes re: data structuring
# 1. I have two sets of columns in the dataset, one that contains zeros when predators 
#-- were absent and one that doesn't contain zeros
# the idea behind this is that if I put zeros in, it will bring the averages way down, and 
#-- affect the diagnostic plots for assumptions. Might rethink this after chatting with Mark

#if I find that there's no difference in scores among treatments, that's probably not great
#want to show that there are equal scores within each treatemnt, but not among treatments
#only issue I forsee are variance differnces among treatments, esp. if certain treatments
#-- can't receive all scores.

#from here, make a df with averages for both sets of intervals

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plyr)

getwd()
#ptl<-read.csv("Data/2019.3.3.ptl.with.interval.factor.csv")
#ptl<-read.csv("Data/2019.3.2.ptl.five.min.intervals.csv")
#ptl<-read.csv("Data/2019.3.2.ptl.ten.min.intervals.csv")
ptl<-read.csv("Data/2019.3.19.ptl.csv")
ptl$Treatment<-ordered(ptl$Treatment, c("Low", "Medium","High"))
ptl$Treatment<-ordered(ptl$Treatment, c("Low", "Medium","High","Control"))
#subsetting for only 60 sec interval data
ptl<-ptl[ptl$Interval=="Ten.minute",]
#subsetting data for trials 4 and 5, and just trial 6
ptl.4.5<-ptl[ptl$Trial<6,]

ptl.6<-ptl[ptl$Trial>5,]
#complete cases df for data where predators were not seen in photo
ptl.4.5.cc<-ptl.4.5[complete.cases(ptl.4.5),]

ptl.6.cc<-ptl.6[complete.cases(ptl.6),]

#NOTE: don't even have to make it into a tibbl, just run raw data
#-------------------#

#but if you wanted to start grouping things, here's how you could do it
ptl.score.int<-ptl %>%
  group_by(Trial,Reef,Treatment,Species,Interval) %>%
  summarize(Score.zero = (mean(Score.zero)))
View(ptl.score.int)

ptl.count.int<-ptl %>%
  group_by(Trial,Reef,Treatment,Species,Interval) %>%
  summarize(Count.zero = (mean(Count.zero)))
View(ptl.count.int)

#five minute intervals:
#averaging scores and counts per treatment, then joining tibbls, using data with zeros
ptl.score.five<-ptl %>%
  group_by(Trial,Reef,Treatment,Species) %>%
  summarize(avg.score = (mean(Score.zero)))
View(ptl.score)

ptl.count.five<-ptl %>%
  group_by(Trial,Reef,Treatment,Species) %>%
  summarize(avg.count = (mean(Count.zero)))
View(ptl.count)

#joining df's
count.and.score.five<-left_join(ptl.score.five,ptl.count.five,by=c("Trial","Reef","Treatment","Species"))
#adding a column to specify the interval read rate (5 min vs. 10 min)
count.and.score.five$Interval<-rep(5,nrow(count.and.score.five))
View(count.and.score.five)

#ten minute interval
ptl.score.ten<-ptl %>%
  group_by(Trial,Reef,Treatment,Species) %>%
  summarize(avg.score = (mean(Score.zero)))
#View(ptl.score)

ptl.count.ten<-ptl %>%
  group_by(Trial,Reef,Treatment,Species) %>%
  summarize(avg.count = (mean(Count.zero)))
#View(ptl.count)

#joining df's

count.and.score.ten<-left_join(ptl.score.ten,ptl.count.ten,by=c("Trial","Reef","Treatment","Species"))
count.and.score.ten$Interval<-rep(10,nrow(count.and.score.five))
View(count.and.score.ten)

#combining two df's
count.and.score.combo<-left_join(count.and.score.five,count.and.score.ten, by=c("Trial","Reef","Treatment","Species","avg.score","avg.count"))
View(count.and.score.combo)

#simple stats with zero's added to data
tapply(ptl$Count.zero,ptl$Treatment,mean)
tapply(ptl$Count.zero,ptl$Treatment,var)

tapply(ptl$Score.zero,ptl$Treatment,mean)
tapply(ptl$Score.zero,ptl$Treatment,var)

#without zeros, using the complete cases only
tapply(ptl.cc$Count,ptl.cc$Treatment,mean)
tapply(ptl.cc$Count,ptl.cc$Treatment,var)

tapply(ptl.cc$Score,ptl.cc$Treatment,mean)
tapply(ptl.cc$Score,ptl.cc$Treatment,var)

#plots####
#counts vs. read numbers to see if 5 mins differ that much from 10 min observations
bargraph.CI(x.factor = Read.number.every.5.min, response = Count.zero, group = Treatment, legend=TRUE, main="Counts over time per treatment", xlab="Read number", ylab="Predator count per reef (mean +/- se)", x.leg=13, yleg=4500, data = ptl)

#using tibbles with average count and score data precalculated
#avg. count and score for 5 min reads
bargraph.CI(x.factor = Treatment, response = avg.count, legend=TRUE, main="5 minute observations, Counts per treatment", xlab="Read number", ylab="Predator count per reef (mean +/- se)", x.leg=13, yleg=4500, data = ptl.count.five)
bargraph.CI(x.factor = Treatment, response = avg.score, group = Treatment, legend=TRUE, main="5 minute observations, scores per treatment", xlab="Read number", ylab="Predator count per reef (mean +/- se)", x.leg=13, yleg=4500, data = ptl.score.five)
#avg. count and score for 10 min reads
bargraph.CI(x.factor = Treatment, response = avg.count, legend=TRUE, main="Counts over time per treatment", xlab="Read number", ylab="Predator count per reef (mean +/- se)", x.leg=13, yleg=4500, data = ptl.count.ten)
bargraph.CI(x.factor = Treatment, response = avg.score, group = Treatment, legend=TRUE, main="Counts over time per treatment", xlab="Read number", ylab="Predator count per reef (mean +/- se)", x.leg=13, yleg=4500, data = ptl.count.ten)

#using raw data and grouping by interval value####
#counts by treatment, I think this is psedoreplicated, so only using the 10 minute interval data here
bargraph.CI(x.factor = Treatment, response = Count.zero, main="Counts per treatment, 10 minute data only", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl)
#counts by intervals
bargraph.CI(x.factor = Treatment, response = Count, group = Interval, legend=TRUE, main="Counts over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)",x.leg = 3, data = ptl)

#no real difference in the number of total fish seen per treatment between two intervals
bargraph.CI(x.factor = Species, response = Count.zero, group = Interval, legend=TRUE, main="Counts over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)", data = ptl)
#no real difference in the number of each species seen between two intervals
#bargraph.CI(x.factor = Treatment, response = Count.zero, group = Species, legend=TRUE, main="Counts over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)", data = ptl)
bargraph.CI(x.factor = Interval, response = Count.zero, group = Species, legend=TRUE, main="Counts over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)", data = ptl)
#looks like the same trends for each species, regardless of 5 vs. 10 minute interval

#scores
#overall by treatment, using only 10-minute interval data
bargraph.CI(x.factor = Treatment, response = Score.zero, main="Scores per treatment, 10 minute data only", xlab="Treatment", ylab="Predator score per reef (mean +/- se)", data = ptl)
#by treatment and interval, comparing two time periods
bargraph.CI(x.factor = Treatment, response = Score.zero, group = Interval, legend=TRUE, main="Scores over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)", data = ptl)
#no real difference in the scores seen per treatment between two intervals
bargraph.CI(x.factor = Species, response = Score.zero, group = Interval, legend=TRUE, main="Scores over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)", data = ptl)
#no real difference in the scores of each species seen between two intervals
#bargraph.CI(x.factor = Treatment, response = Score.zero, group = Species, legend=TRUE, main="Scores over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)", data = ptl)
bargraph.CI(x.factor = Interval, response = Score.zero, group = Species, legend=TRUE, main="Scores over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)", data = ptl)
#looks like the same trends for each species, regardless of 5 vs. 10 minute interval

#linear modeling##############
#I want to be able to compare 10-min intervals counts/scores
#--to 5-mon interval counts/scores to see if there is a sig. diff.

mod1<-lm(Count~Treatment*Interval, data=ptl)
hist(resid(mod1))
qqnorm(resid(mod1))
#boxplot(Score~Treatment, data=ptl)
anova(mod1)
TukeyHSD(mod1,ptl$Count~ptl$Treatment,ordered=TRUE)

mod1<-aov(Count~Treatment, data=ptl)
hist(resid(mod1))
qqnorm(resid(mod1))
#boxplot(Score~Treatment, data=ptl)
anova(mod1)
TukeyHSD(mod1,ordered=TRUE)

mod2<-lm(Score.zero~Treatment*Interval,data=ptl)
hist(resid(mod2))
qqnorm(resid(mod2))
#boxplot(Score~Treatment, data=ptl)
anova(mod1)

#boxplots for quick comparison

#2018 t4.5

#using data with zeros:
boxplot(ptl.4.5$Count.zero~ptl.4.5$Treatment)
boxplot(ptl.4.5$Score.zero~ptl.4.5$Treatment)

bargraph.CI(x.factor = Treatment, response = Count.zero, main="T4 + 5 Count data, with zeros", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl.4.5)
bargraph.CI(x.factor = Treatment, response = Score.zero, main="T4 + 5 Scores data, with zeros", xlab="Treatment", ylab="Predator score per reef (mean +/- se)", data = ptl.4.5)

#using without zeros, and only complete cases (t4 and 5)
boxplot(ptl.4.5.cc$Count~ptl.4.5.cc$Treatment)
boxplot(ptl.4.5.cc$Score~ptl.4.5.cc$Treatment)

bargraph.CI(x.factor = Treatment, response = Count, main="T4 + 5 Count data, no zeros", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl.4.5.cc)
bargraph.CI(x.factor = Treatment, response = Score, main="T4 + 5 Score data, no zeros", xlab="Treatment", ylab="Predator score per reef (mean +/- se)", data = ptl.4.5.cc)

#trial 6 data only
boxplot(ptl.6$Count.zero~ptl.6$Treatment)
boxplot(ptl.6$Score.zero~ptl.6$Treatment)

bargraph.CI(x.factor = Treatment, response = Count.zero, main="T6 Count data, with zeros", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl.6)
bargraph.CI(x.factor = Treatment, response = Score.zero, main="T6 Scores data, with zeros", xlab="Treatment", ylab="Predator score per reef (mean +/- se)", data = ptl.6)

#using without zeros, and only complete cases (t4 and 5)
boxplot(ptl.6.cc$Count~ptl.6.cc$Treatment)
boxplot(ptl.6.cc$Score~ptl.6.cc$Treatment)

bargraph.CI(x.factor = Treatment, response = Count, main="T6 Count data, no zeros", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl.6.cc)
bargraph.CI(x.factor = Treatment, response = Score, main="T6 Score data, no zeros", xlab="Treatment", ylab="Predator score per reef (mean +/- se)", data = ptl.6.cc)

#2019.3.18, not quite sure what the figs indicate at this point, might have to pick up slack 
#-- on other treatments for trial 6, because it'stough to separate trials out with only
#-- high and control treatments being analyzed

#on one hand, it shows that high and control didn't have any differences in predator activity,
#-- at least not when zeros were excluded

#revisualizing data with binary plots to see the proportion of data that are zero

# proportion of data that are ==, >, <, a certain value
plot(factor(Score.zero==0)~Treatment, data=ptl.4.5)
plot(factor(Count.zero >1)~Treatment, data=ptl.4.5)
#shows that there are a ton of 0's in the data, more than 90% of all the data are 0
#not usable?

#binomial analyses for zero-inflated data:

#code from: https://stats.stackexchange.com/questions/156643/methods-to-analyze-zero-inflated-data

#stock code:

#d2 <- subset(d, type != "B")
#d2$type <- factor(d2$type)
#m1 <- glm(factor(intact > 0) ~ type, data = d2,
          #family = binomial(link = "probit"))
#summary(m1)

#trials 4 and 5
d2<-subset(ptl,Trial<6) #only looking at t4.5 data

#scores from t4.5
m1 <- glm(factor(Score.zero) ~ Treatment, data = d2,
          family = binomial(link = "probit"))
summary(m1)
Anova(m1,type="II")
hist(resid(m1))
qqnorm(resid(m1))

#counts from t4.5

m2 <- glm(factor(Count.zero) ~ Treatment, data = d2,
            family = binomial(link = "probit"))
summary(m2)
Anova(m2,type="II")
hist(resid(m2))
qqnorm(resid(m2))

#permanova (1st attempt on 2019.3.19)
library(vegan)

## Using plyr to reorder df for permanova
arrange(score, sex,y)
ddply(score, c('sex', 'y'))

d2<-arrange(d2, Treatment)
#ddply(ptl, c('sex', 'y'))
View(d2)

#Univariate permanova
#extracting the columns with invert data only 

qqnorm(invert$purple.urchin)
qqline(invert$purple.urchin)

#two-way PERMANOVA MODEL observing site,protection and the interaction with euclidean measure because it won't run with Bray Curtis( dosen't like the zeroes). 

permanovamodel2<-adonis(purple.urchin~MPA*region, data = invert, permutations = 999, method="euclidean", by= "terms")
permanovamodel2

permanovamodel2<-adonis(Count.zero~Treatment, data = d2, permutations = 999, method="euclidean", by= "terms")
permanovamodel2

###permdisp####
# test of homogeneity of variances (PERMDISP - multivariate analogue of Levene's Test)
mv_dist<-vegdist(invert$purple.urchin, method = "euclidean")
mv_dist
# define location groups to compare homogeneity of dispersion
groups2<-factor(c(rep(1,36),rep(2,36)),labels=c("Yes","No"))

# Calculate multivariate dispersions
mod<-betadisper(mv_dist,groups2)
mod
anova(mod)

#code from E. nava

adonis(Count.zero~Treatment, data=d2, permutations=999)

#I'm not sure this is working, because I'm getting the same answer for
#-- scores and count data

#explanation(I think)
#What this is doing is telling me if there's a sig difference in the number of
#-- zeros by treatment


#trial 6
d3<-subset(ptl,Trial>5)

#scores from t6
m3 <- glm(factor(Score.zero>0) ~ Treatment, data = d3,
            family = binomial(link = "probit"))
summary(m3)
Anova(m3,type="II")
hist(resid(m3))
qqnorm(resid(m3))

#counts from t6
m4 <- glm(factor(Count.zero > 0) ~ Treatment, data = d3,
          family = binomial(link = "probit"))
summary(m4)
Anova(m4,type="II")
hist(resid(m4))
qqnorm(resid(m4))

