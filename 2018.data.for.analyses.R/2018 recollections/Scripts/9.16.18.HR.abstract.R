##############################################################################
# recollections data revamped with one row for each reef                     #
#   George Jarvis                                                            #
#    9/16/2018                                                               #
##############################################################################

rm(list=ls())

#load packages
library(sciplot)
library(lme4)
library(lmerTest)
library(car)

getwd()
sc<-read.csv("Data/9.16.18.recollections.data.csv")
sc$trial<-as.factor(sc$trial)#made trial a factor
#only including adults that we tagged initially (immatures were listed as NA in datasheet)
sc<- sc[complete.cases(sc), ]#25 rows of data
#sc$Growth<-as.numeric(sc$Growth)
head(sc)

sc.mod.r<-lm(prop.male~treatment, data=sc)
hist(resid(sc.mod.r))
qqnorm(resid(sc.mod.r))
boxplot(resid(sc.mod.r))
boxplot(prop.male~treatment, data=sc)

anova(sc.mod.r)
summary(sc.mod.r)

#pooling all trials together

#plotting proportion recollected by risk
#1) males
bargraph.CI(x.factor = treatment, response = prop.male, main="proportion male", data = sc)
#seems like no sig difference in growth among treatments
#2)females
bargraph.CI(x.factor = treatment, response = prop.female, main="proportion female", data = sc)
#seems like no sig difference in growth among treatments
#3)now want to look at male:female ratio/treatment
  #only want complete cases, which are ones where there was at least 1 male and 1 female recollected per reef
  #otherwise, there's no way to get a single ratio for each reef

sc<- sc[complete.cases(sc), ]#25 rows of data

sc.mod.sr<-lm(final.sex.ratio~treatment, data=sc)
hist(resid(sc.mod.sr))
qqnorm(resid(sc.mod.sr))
boxplot(resid(sc.mod.sr))
boxplot(final.sex.ratio~treatment, data=sc)
anova(sc.mod.sr)
summary(sc.mod.sr)
#plotting
bargraph.CI(x.factor = treatment, response = final.sex.ratio, main="final sex ratio male:fem", data = sc)
#will revisit this, but for now, it seems like there were very high numbers of males
#recollected from low and med risk treatments, and fewer in the low-risk trt.
#these recollection ratios are way higher than the natural ratios (1:3, m:f)

#there must be a way to normalize this by the number of total fish that were recollected
#going to try and divide these final sex ratios by the total number recollected

sc.mod.srn<-lm(norm.sex.ratio~treatment, data=sc)
hist(resid(sc.mod.srn))
qqnorm(resid(sc.mod.srn))
boxplot(resid(sc.mod.srn))
boxplot(norm.sex.ratio~treatment, data=sc)
anova(sc.mod.srn)
summary(sc.mod.srn)
bargraph.CI(x.factor = treatment, response = norm.sex.ratio, main="normalized sex ratio male:fem", data = sc)

#also want to run this again as a model that includes the nubmer of fish recollected
#as a covariate 

sc.mod.srt<-lm(final.sex.ratio~treatment+total, data=sc)
hist(resid(sc.mod.srt))
qqnorm(resid(sc.mod.srt))
boxplot(resid(sc.mod.srt))
boxplot(final.sex.ratio~treatment, data=sc)
anova(sc.mod.srt)
summary(sc.mod.srt)
#not better....I might have to go through a whole model selection process (non-linear)