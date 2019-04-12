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
library(vegan)
library(psych)
library(MASS)
library(boot)
library(effects)
library(simpleboot)

getwd()
#ptl<-read.csv("Data/2019.3.3.ptl.with.interval.factor.csv")
#ptl<-read.csv("Data/2019.3.2.ptl.five.min.intervals.csv")
#ptl<-read.csv("Data/2019.3.2.ptl.ten.min.intervals.csv")
ptl<-read.csv("Data/2019.4.4.ptl.csv")#these data don't include the super rare species
#   I got reid of CAPR, HYRU, PANE, and GINI
ptl<-ptl[,1:15]
View(ptl)
ptl$Treatment<-ordered(ptl$Treatment, c("Low", "Medium","High"))
ptl$Treatment<-ordered(ptl$Treatment, c("Low", "Medium","High","Control"))
#subsetting for only 60 sec interval data
ptl<-ptl[ptl$Interval=="Ten.minute",]
#subsetting data for trials 4 and 5, and just trial 6
ptl<-ptl[ptl$Trial<6,]
View(ptl.4.5)

ptl.6<-ptl[ptl$Trial>5,]
#complete cases df for data where predators were not seen in photo
ptl.4.5.cc<-ptl.4.5[complete.cases(ptl.4.5),]
ptl<-subset(ptl,Treatment!="Control")

#removing control treatment from df

ptl.6.cc<-ptl.6[complete.cases(ptl.6),]

#NOTE: don't even have to make it into a tibbl, just run raw data
#-------------------#

#2019.3.20

ptl.score.int.a<-ptl %>%
  na.omit() %>%
  group_by(Trial,Reef,Treatment,Species,Interval) %>%
  summarize(Count.zero = (mean(Count.zero)))
View(ptl.score.int.a)

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
boxplot(ptl.4.5.cc$Count~ptl.4.5.cc$Treatment, xlab="Treatment", ylab="Count", main="Trial 4 and 5 predator abundance, when predators captured in frame")
#I think the plot for counts works well as a barplot, not as a boxplot
boxplot(ptl.4.5.cc$Score~ptl.4.5.cc$Treatment, xlab="Treatment", ylab="Score", main="Trial 4 and 5 predator scores, when predators present")


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
plot(factor(Count.zero>1)~Treatment, data=ptl.4.5)
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

ptl<-arrange(ptl, Treatment)
#ddply(ptl, c('sex', 'y'))
View(ptl)

d2<-arrange(d2,Treatment)

#Univariate permanova
#extracting the columns with invert data only 

qqnorm(invert$purple.urchin)
qqline(invert$purple.urchin)

#two-way PERMANOVA MODEL observing site,protection and the interaction with euclidean measure because it won't run with Bray Curtis( dosen't like the zeroes). 

#permanovamodel2<-adonis(purple.urchin~MPA*region, data = invert, permutations = 999, method="euclidean", by= "terms")
#permanovamodel2

#getting rid of NA's

#tells me how many na there is within each column
colSums(is.na(d2))
#omit na within score column
d2<-d2[!is.na(d2$Score.zero),]

#removing control treatment from df
d2<-subset(ptl,Treatment!="Control")

permanovamodel2<-adonis(Score.zero~Treatment, data = d2, permutations = 999, method="euclidean", by= "terms")
permanovamodel2

#figuring out how many rows of each treatment there are

#need the sample size for each level of treatment
count(d2$Treatment)

library(psych)

###permdisp####
# test of homogeneity of variances (PERMDISP - multivariate analogue of Levene's Test)
#mv_dist<-vegdist(invert$purple.urchin, method = "euclidean")
#mv_dist

mv_disp<-vegdist(d2$Score.zero, method = "euclidean")
mv_disp
# define location groups to compare homogeneity of dispersion
groups2<-factor(c(rep(1,950),rep(2,841),rep(3,350)),labels=c("High","Low","Medium"))

# Calculate multivariate dispersions
mod<-betadisper(mv_disp,groups2)
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

#griffin's code
IND_data<-read.csv("ind_survey_V2.csv")

#### Ordinal logistic regression ####
IND_ordinal <- polr(factor(richness) ~ Height, data = IND_data) # model 
summary(IND_ordinal) ## Summary, but with no p-value, if you need a p-value you will have to compare t-value against normal distribution. See below....

ctable <- coef(summary(IND_ordinal)) ## put the coefficients in a table
p_IND_ordinal<- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2  ## calcuate p-values
ctable <- cbind(ctable, "p value" = p_IND_ordinal) ## p-values to table
ctable # presto, p-values

#my version, not sure this is what I want

ptl.ordinal <- polr(factor(Score) ~ Treatment, data = ptl.4.5.cc) # model 
summary(IND_ordinal) ## Summary, but with no p-value, if you need a p-value you will have to compare t-value against normal distribution. See below....

ctable <- coef(summary(ptl.ordinal)) ## put the coefficients in a table
p_IND_ordinal<- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2  ## calcuate p-values
ctable <- cbind(ctable, "p value" = p_IND_ordinal) ## p-values to table
ctable # presto, p-values

#restructuring data to show when predators were present,
# the proportion of time that they were within distance of being perceived
# as a threat by gobies (score >3)

plot(factor(Score.zero==0)~Treatment, data=ptl) #about 20% of the data have scores equal to zero
plot(factor(Score.zero>3)~Treatment, data=ptl) #about 8% of fish seen in HR treatemnt were considered a threat
                                              #less than 5% of fish in med risk were seen as a threat

plot(factor(Score==0)~Treatment, data=ptl)
plot(factor(Count.zero>3)~Treatment, data=ptl)

#making new variable for when predators were close enough to be perceived
# as risky by gobies

library(dplyr)

ptl.4.5.new<-ptl.4.5.cc %>%
  mutate(newvar = ifelse(factor1 > 5 & 
                           factor2 < 19 & 
                            & 
                           !is.na(factor4), 1, NA))

library(dplyr)

dataset %>%
  mutate(newvar = ifelse(factor1 > 5 & 
                           factor2 < 19 & 
                           (factor3=="b" | factor3=="c") & 
                           !is.na(factor4), 1, NA))

#try this
egg.den.bio <- egg.den.bio %>% 
  mutate(Year=ifelse((Trial<4),2017,2018))

#trying with ptl data (including zeros)
#score of >3 = threat, score <4 = no threat
ptl.4.5.new <- ptl %>% 
  mutate(prop.threat=ifelse((Score.zero>3),1,0))
ptl.4.5.new <- ptl %>% 
  mutate(prop.threat=ifelse((Score>3),1,0))

View(ptl.4.5.new)

bargraph.CI(x.factor = Treatment, response = prop.threat, main="T4 + 5 Count data, with zeros", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl.4.5.new)
bargraph.CI(x.factor = Treatment, response = Score.zero, main="T4 + 5 Scores data, with zeros", xlab="Treatment", ylab="Predator score per reef (mean +/- se)", data = ptl.4.5)

ptl.4.5.cc <- ptl.4.5.cc %>% 
  mutate(prop.threat=ifelse((Score>3) &
                              !is.na(Score),1,0))
View(ptl.4.5.cc)


#quick plotting for species
bargraph.CI(x.factor = Treatment, response = Score.zero, group = Species, legend=TRUE, main="Scores over time per treatment and time interval, raw data", xlab="Read number", ylab="Predator count per reef (mean +/- se)", data = ptl)
tapply(ptl$Count.zero,ptl$Species,mean,na.rm=TRUE)
#plotting with new df that excludes rare species
bargraph.CI(x.factor = Treatment, response = prop.threat, main="T4 + 5 Count data, with zeros", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl.4.5.new)
bargraph.CI(x.factor = Treatment, response = Score, main="T4 + 5 Scores data, with zeros", xlab="Treatment", ylab="Predator score per reef (mean +/- se)", data = ptl)

mod1<-glm(cbind(presence, absence) ~ 1 + treatment + year, family=binomial)

#removed low treatment because not achievable for predator to be berceived as risk here
ptl.no.low<-subset(ptl.4.5.new,Treatment!="Low")

mod1<-glm(cbind(prop.threat, absence) ~ 1 + treatment + year, family=binomial)

mod1<-glm(prop.threat~Treatment,family=binomial,data=ptl.no.low)
qqnorm(resid(mod1))
hist(resid(mod1))
summary(mod1)
Anova(mod1)

#proportion of photos that don't contain predators
bargraph.CI(x.factor = Treatment, response = Count, main="T4 + 5 Count data, with zeros", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl)
bargraph.CI(x.factor = Treatment, response = Count.zero, main="T4 + 5 Count data, with zeros", xlab="Treatment", ylab="Predator count per reef (mean +/- se)", data = ptl)

mod1<-lm(Count.zero~Treatment,family=poisson,data=ptl)
qqnorm(resid(mod1))
hist(resid(mod1))

#plotting without low treatment
#proportion of predators perceived as a threat
#of all photos taken that contained predators
bargraph.CI(x.factor = Treatment, response = (ptl$Score>3), main="prop perceived", xlab="Treatment", ylab="Proportion of predators perceived as threat", data = ptl)
#of all photos, including zeros
bargraph.CI(x.factor = Treatment, response = prop.threat, main="prop perceived", xlab="Treatment", ylab="Proportion of predators perceived as threat", data = ptl)

#prop of photos that captured predators at all
bargraph.CI(x.factor = Treatment, response = count.prop, main="prop photos without predators", xlab="Treatment", ylab="Proportion of photos with predators", data = ptl)

modc<-glm(count.prop~Treatment,family=poisson,data=ptl)
summary(modc)
Anova(modc, type="III")
hist(resid(modc))
qqnorm(resid(modc))

modd<-glm(prop.threat~Treatment,family=binomial,data=ptl)
Anova(modd, type="III")
hist(resid(modd))

lmodel <- lm(Ozone ~ Wind)
lboot <- lm.boot(lmodel, R = 1000)
summary(lboot)

count.boot<-lm.boot(modc,R=1000)

#bootstrapping 95% confidence intervals with boot package

m1 <- lm(Fertility ~ ., swiss)
betahat.boot <- Boot(mod1, R=999) # 199 bootstrap samples--too small to be useful
summary(betahat.boot)  # default summary
confint(betahat.boot)
hist(betahat.boot)
boot.ci(betahat.boot, type=c("stud","norm"))

boot.ci(mod1, type=c("stud","norm"))

m1 <- lm(Fertility ~ ., swiss)
betahat.boot <- Boot(m1, R=199) # 199 bootstrap samples--too small to be useful
summary(betahat.boot)  # default summary
confint(betahat.boot)
hist(betahat.boot)

#manually making confidence intervals based on calcualted means of data
#calculating means


a <- 5
s <- 2
n <- 20
error <- qnorm(0.975)*s/sqrt(n)
left <- a-error
right <- a+error
left
#[1] 4.123477
right
#[1] 5.876523

#kathryn's code
View(stacked_groups)
combined_algae<-data.frame(cbind(SAHO,SAPA,DIUN,ZOFA))
summary(combined_algae)
stacked_groups<-stack(combined_algae)
stacked_groups
invert.mod<-aov(log(values)~ind, data = stacked_groups)
summary(invert.mod)
TukeyHSD(invert.mod)

qqnorm(resid(invert.mod))
qqline(resid(invert.mod))
hist(resid(invert.mod))
?effect
k<-effect

ef<-effect("Treatment", modc)
ef
df<-data.frame(ef)
df
df$Treatment<-ordered(df$Treatment, c("Low", "Medium","High"))
x$fit<-exp(x$fit)
x$lower<-exp(x$lower)
x$upper<-exp(x$upper)
x$se<-exp(x$se)
x

#This is how to change the names of the x axis bars and to italicize them
mylabels<-c(expression(paste(italic("D. undulata"))),
            expression(paste(italic("S. horneri"))),
            expression(paste(italic("S. palmeri"))),
            expression(paste(italic("Z. farlowii"))))


###graph for inverts on each algal species
theme_update(plot.title = element_text(hjust = 0.5)) #this centers the title 
ggplot(data = x, aes(x=ind, y=fit, fill = ind)) + #your x and y, x=algal species, y = means of forgaing rate for each species
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  scale_x_discrete(labels=mylabels)+ #This scale_x_discrete part will change the x axis bar names for the labels and the limits changes the order of the bars
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower,ymax = upper), position = position_dodge(.9), 
                width = 0.2) + coord_cartesian (ylim = c(0,1)) + ### grapging your standard errors using confidence intervals 
  ggtitle('Invertebrate Abundance per Alga') + # title of your graph
  theme(plot.title = element_text(size=18, face = "bold", family = "sans"))+ # changing the plot title size and font
  xlab("Algal Species") + #title for the x axis
  ylab(expression(paste("Invertebrate Abundance"))) + #title for the y axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ #takes away the gray grid in the background 
  theme(axis.title.x = element_text(size=18, family = "sans"),#this changes the size and angle of the x-axis column labels
        axis.text.x  = element_text(angle=0, vjust=0.5, size=12)) + #this changes size and font of x axis
  theme(axis.title.y = element_text(size = 18, family = "sans"),#this changes the size and angle of the x-axis column labels
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12))+ #this changes size and font of y axis
  scale_fill_manual(values = pal(4))+
  theme(legend.position = "none") #this deletes the legend
#theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ##this takes away the axis ticks, personal preference

plot(invert.mod)

###plotting###

#counts, but with reef as replicate
count.means<-with(ptl, aggregate((count.prop), list(Treatment=Treatment,Reef=Reef,Species=Species), mean))
#now apply the se function to the 4th column [,3]
count.means$se<-with(ptl, aggregate((count.prop), list(Treatment=Treatment,Reef=Reef,Species=Species), function(x) sd(x)/sqrt(length(x))))[,4]
count.means

#down the rabbit hole we go...
count.means1<-with(count.means, aggregate((x), list(Treatment=Treatment), mean))
count.means1$se<-with(ptl, aggregate((count.prop), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
count.means1

mod1<-lm(x~Treatment,data=count.means)
hist(resid(mod1))
qqnorm(resid(mod1))
summary(mod1)
Anova(mod1,type="III")

png(filename = "Output/counts.treatment.new.png", width = 900, height = 800)

c<- ggplot(count.means1, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=25, colour="black"), axis.title=element_text(size=30,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 35, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
c + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                            position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Photos with Predators")

dev.off()

#scores

mod2<-aov(x~Treatment,data=score.means)
hist(resid(mod2))
qqnorm(resid(mod2))
summary(mod2)
Anova(mod2)
TukeyHSD(mod2,score.means$Treatment)
TukeyHSD(mod2, "Treatment", ordered = TRUE)

#scores, but with reef as replicate
score.means<-with(ptl, aggregate((prop.threat), list(Treatment=Treatment,Reef=Reef,Species=Species), mean))
#now apply the se function to the 4th column [,3]
score.means$se<-with(ptl, aggregate((prop.threat), list(Treatment=Treatment,Reef=Reef,Species=Species), function(x) sd(x)/sqrt(length(x))))[,4]
score.means

#down the rabbit hole we go...
score.means1<-with(score.means, aggregate((x), list(Treatment=Treatment), mean))
score.means1$se<-with(score.means, aggregate((x), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
score.means1
score.means1Treatment<-ordered(score.means1$Treatment, c("Low", "Medium","High"))

mod1<-aov(x~Treatment,data=score.means)
hist(resid(mod1))
qqnorm(resid(mod1))
summary(mod1)
Anova(mod1,type= "II")
TukeyHSD(mod1, "Treatment", ordered = TRUE)

png(filename = "Output/scores.treatment.png", width = 900, height = 800)

c<- ggplot(score.means1, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=25, colour="black"),axis.text.y=element_text(size=25, colour="black"), axis.title=element_text(size=30,face="bold")) +
  theme(axis.title.y = element_text(size= 30, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
c + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                   position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Predators Perceived as Threat")

dev.off()


###trying again, but now added function for CI
#scores, but with reef as replicate
score.means<-with(ptl, aggregate((prop.threat), list(Treatment=Treatment,Reef=Reef,Species=Species), mean))
#now apply the se function to the 4th column [,3]
score.means$CI<-with(ptl, aggregate((prop.threat), list(Treatment=Treatment,Reef=Reef,Species=Species), function(x) qnorm(0.975)*sd(x)/sqrt(length(x))))[,4]
score.means

#down the rabbit hole we go...
score.means1<-with(score.means, aggregate((x), list(Treatment=Treatment), mean))
score.means1$se<-with(score.means, aggregate((x), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
score.means1
score.means1Treatment<-ordered(score.means1$Treatment, c("Low", "Medium","High"))

mod1<-aov(x~Treatment,data=score.means)
hist(resid(mod1))
qqnorm(resid(mod1))
summary(mod1)
Anova(mod1,type= "II")
TukeyHSD(mod1, "Treatment", ordered = TRUE)

png(filename = "Output/scores.treatment.png", width = 900, height = 800)

c<- ggplot(score.means1, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=25, colour="black"),axis.text.y=element_text(size=25, colour="black"), axis.title=element_text(size=30,face="bold")) +
  theme(axis.title.y = element_text(size= 30, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
c + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                   position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Predators Perceived as Threat")

dev.off()

















#trying again

count.means$upper<-df$upper
count.means$lower<-df$lower
count.means

count<-ptl %>%
  group_by(Species) %>%
  summarize(mean = (mean(count.prop)))
View(count)


View(ptl)
