# Description: Analyses and plots for predator behavior in Jarvis and Steele
# Author: George C Jarvis
# Date: Sun May 31 10:38:29 2020
# Notes: Re: coding for predator activity:
# 1) Time-lapses for each reef were subsampled, and 10-20 photos were analyzed per time lapse, depending on its duration
# 2) Each photo received a single score for predator presence (0 if no predators present, 1 if any predators present)
# 3) If predators were present, they were scored based on their location, relative to reef structure
#   a. if far from structure, it was marked as present, but not as either lethal or sublethal (i.e. low threat to gobies)
#   b. if close to structure but not close enough to consume gobies (i.e. in medium risk, high-risk caged, or high-risk uncaged
#   - treatments only), then the predator was marked as a 1 for sublethal threat
#   c. if close to structure and close enough to consume gobies (i.e. in high-risk caged and high-risk uncaged
#   - treatments only), then the predator was marked as a 1 for lethal threat
# 4) If there were multiple predators in the photo, and they were seen in different locations relative to the reef, then 
#   - each category for predator scores could have received a value of 1
# 5) Each category (presence, sublethal, lethal) received a score of either 1 or 0 for each photo
# --------------

rm(list=ls())

#packages####
library(lme4)
library(lmerTest)
library(ggplot2)
library(tidyverse)

#type III ANOVA
options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#importing data####
pred<-read.csv(file = "Data/2020_5_31_predator_raw_data.csv")
View(pred)
#removing NA's
pred.rm.na<-na.omit(pred)
head(pred.rm.na)

#rename "Treatment.combo" to "Treatment"
pred.rm.na<-rename(pred.rm.na, Treatment = Treatment.combo)
pred.rm.na<-rename(pred.rm.na, T6.comparison = Treatment.t6)

#adding two columns for year, one for year (int), and one for year.fact (fact) to df;
# - as a proxy for tagging procedure, where trials 1-3 = 2017, and 4-6 = 2018
pred.rm.na$Year <- ifelse(pred.rm.na$Trial <=3, 2017, 2018)

#making year and trial factors for analysis
pred.rm.na$Year<- as.factor(pred.rm.na$Year)
pred.rm.na$Trial<-as.factor(pred.rm.na$Trial)

# exporting data in wide format
# write.csv(pred.rm.na,"Data\\2019.1.2.predator.raw.data.csv", row.names = FALSE)

# importing egg count data for avg. number of gobies inhabiting each reef
# Will add this to predator data when I wrangle them into wide format
repro<-read.csv(file= "Data/egg.counts.2019.12.23.csv")
#making Year.fact a factor
repro$Year<-as.factor(repro$Year)
repro$Trial<-as.factor(repro$Trial)


#data wrangling####

#wide format (for analyses)####

# want the proportion of photos per reef +/- predators, 
# - predators close enough to be a sublethal threat, 
# - and predators close enough to be a lethal threat

# Grouping by Year, Trial, Reef, Treatment, and Treatment.t6 (for t6 comparison) 
# Read number is the level of observation
# Each reef is a replicate, so that is what will be used to calculate means and SEM for plots

p<-with(pred.rm.na, aggregate(list(contained.pred,sublethal.threat,lethal.threat), 
                              list(Trial=Trial,Reef=Reef,Treatment=Treatment,
                                   T6.comparison=T6.comparison,Year=Year), mean))
#View(p)

#Column names for columsn 6, 7 and 8 are not correct, need to change them to "Present", "Sublethal.Threat", 
# and "Lethal.Threat", respectively

p <- p %>% 
  rename(Present = 6, Sublethal.Threat = 7, Lethal.Threat = 8)
#View(p)

# Bringing in avg. numer of gobies inhabiting each reef as well, adding it to "p" df

pr<-left_join(p,repro,by=c("Trial","Reef","Treatment","T6.comparison","Year"))
View(pr)

#want to drop unnecessary columns from repro df, i.e. everything but avg.inhab
# they are columns 8:12 and 14:16
pr1<-pr[-c(10:14,16)]
View(pr1)

#long format for plotting#### 
# p-long = "pl"

#for all trials
pr1<-as_tibble(pr1)
pr1

pl<-pr1 %>% gather(Predator.class, Score, Present:Lethal.Threat)
View(pl)

# Creating separate df named "HR.comp" to compare responses between HR-caged and -uncaged treatments in Trial 6
# Will not include year or trial in the model for analyses
# I only did time lapses on high-risk caged and uncaged treatments for T6

HR.comp<-p[p$Trial==6,]
View(HR.comp)

#removing Trial and Year columns

HR.comp<-HR.comp %>% select (-c(Trial,Year))

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

#2020.1.30.edit, just comparing HR caged and uncaged plots with a t-test for each distinction
#1) present
#2) sublethal threat
#3) lethal threat

#loading data (in proper format for t-test)
hr.t.test<-read.csv(file = "Data/2020.1.30.ptl.t.test.cage.artifacts.csv") #sending this to Mark
View(hr.t.test)
#1) present
t.test(hr.t.test$present.high,hr.t.test$present.u)
#no diff in number of predators present

#2) sublethal
t.test(hr.t.test$sublethal.h,hr.t.test$sublethal.u)
#no diff in prop of sublethal

#3) lethal
t.test(hr.t.test$lethal.h,hr.t.test$lethal.u)
#no diff in lethal, in fact, equal means for each treatment

# I think there are some issues here with sample size, because
# - there is so much error and such litle sample sizes for each (n=4)
# - that it would be hard to detect differences, but based on these
# - results it suggests that there were no effect of caging artifacts,
# - at least in the high-risk treatments

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

#ggplot plotting with data in long format###
#first working with data from trials 1-5: proportions of photos where predators were
# - present (low perceived threat, no actual threat), likely perceived as a sublethal threat
# - (high perceived threat, no actual threat), and likely perceived as a high threat (high perceived
# - threat, actual threat)

View(trial.1.5.long)

#formatting data to work with in ggplot

#calculating means and se for each category for proximity
# want a df that contains average score (proxy for proportion of photos that met the criteria)
# - grouped by treatment and predator class

pred1.5<-with(trial.1.5.long, aggregate((Score),list(Treatment=Treatment,Predator.class=Predator.class),mean))
#View(pred1.5)
pred1.5$se<-with(trial.1.5.long, aggregate((Score),list(Treatment=Treatment,Predator.class=Predator.class), 
                                           function(x) sd(x)/sqrt(length(x))))[,3]
#ordering Predator.class values
pred1.5$Predator.class<-ordered(pred1.5$Predator.class, c("Present","Sublethal.Threat","Lethal.Threat"))

#View(pred1.5)

#ggplot master code

#this is the plot that I used to make the black and grayscale figure, with no figure legend, for the MS

png("Output/2019.2.6.9.5x5.5.300dpi.png", width = 9.5, height = 5.5, units = 'in', res = 300)

t1.5.plot<- ggplot(pred1.5, aes(x=Treatment, y=x, fill=Predator.class)) +
  geom_bar(stat="identity", colour= "black", width = 0.9, position="dodge", show.legend = FALSE)+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y = "Proportion of Photos")+ 
  scale_fill_manual(values=c("#666666", "#999999", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.605))
t1.5.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.90)) + theme(text = element_text(family="Arial"))
dev.off()

scale_fill_discrete(name="Predator Classification", labels=c("Present, low threat", "High perceived, no actual threat", 
                                                             "High perceived and actual threat"), values=c("black", "#666666", "grey"))
#messing around with order of fill, have to do this in two steps, first make black and white with no legend, then color (fill discrete)
# - with labels that are correct
# Then have to bring into PPT and change the color of the boxes in the legend to match the colors in the B&W figure


#this is the plot that I used to make a color figure, with a figure legend (see scale fill discrete), for the MS
png("Output/2019.2.6.9.5x5.5.color.300dpi.png", width = 9.5, height = 5.5, units = 'in', res = 300)

t1.5.plot<- ggplot(pred1.5, aes(x=Treatment, y=x, fill=Predator.class)) +
  geom_bar(stat="identity", colour= "black", width = 0.9, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y = "Proportion of Photos")+ 
  theme(legend.position="right") + 
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.605))+
  scale_fill_discrete(name="Threat from Predators", labels=c("Low perceived, no actual threat", "High perceived, no actual threat", 
                                                             "High perceived, actual threat"))
t1.5.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.90)) + theme(text = element_text(family="Arial"))

dev.off()

#now have to do the same thing, but with Trial 6 data only (HR.long df)
#using "T6.comparison" for treatment factor this time

#View(HR.long)

pred.6<-with(HR.long, aggregate((Score),list(T6.comparison=T6.comparison,Predator.class=Predator.class),mean))
#View(pred.6)
pred.6$se<-with(HR.long, aggregate((Score),list(T6.comparison=T6.comparison,Predator.class=Predator.class), 
                                   function(x) sd(x)/sqrt(length(x))))[,3]
#ordering Predator.class values
pred.6$Predator.class<-ordered(pred.6$Predator.class, c("Present","Sublethal.Threat","Lethal.Threat"))

#renaming the labels for high-risk treatment t6.comparison with base R from High to "High - Caged" and "Uncaged" to "High - Uncaged" 
#will eventually have a second figure that has "Risk Treatment" as the x axis, and can just compare caging effects, that will be more clear anyway
#doing it in the new df ("pred.6"), not in the original df used to make calculations ("HR.long")

# Rename by name: change "High" to "Caged"
levels(pred.6$T6.comparison)[levels(pred.6$T6.comparison)=="High"] <- "High - Caged"
levels(pred.6$T6.comparison)[levels(pred.6$T6.comparison)=="Uncaged"] <- "High - Uncaged"

#note: this is the code that I used to make the figure for the MS (see PPT for figure manipulations)

png("Output/2019.2.6.9.5x5.5.trial6.300dpi.png", width = 9.5, height = 5.5, units = 'in', res = 300)

t1.5.plot<- ggplot(pred.6, aes(x=T6.comparison, y=x, fill=Predator.class)) +
  geom_bar(stat="identity", colour= "black", width = 0.60, position="dodge", show.legend = FALSE)+ 
  scale_x_discrete(limits=c("High - Caged", "High - Uncaged"))+
  theme_classic() + 
  labs(x="Risk Treatment", y = "Proportion of Photos")+ 
  scale_fill_manual(values=c("#666666", "#999999", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.82))
t1.5.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.60)) + theme(text = element_text(family="Arial"))
dev.off()
