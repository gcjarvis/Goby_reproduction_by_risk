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
library(FSA) #summarizing samples and means
library(emmeans)

#importing data####
pred<-read.csv(file = "Data/2020_5_31_predator_raw_data.csv")
#View(pred)
#removing NA's
pred.rm.na<-na.omit(pred)
head(pred.rm.na)
View(pred.rm.na)

#rename "Treatment.combo" to "Treatment"
pred.rm.na<-rename(pred.rm.na, Treatment = Treatment.combo)
pred.rm.na<-rename(pred.rm.na, T6.comparison = Treatment.t6)

#adding two columns for year, one for year (int), and one for year.fact (fact) to df;
# - as a proxy for tagging procedure, where trials 1-3 = 2017, and 4-6 = 2018
pred.rm.na$Year <- ifelse(pred.rm.na$Trial <=3, 2017, 2018)

#making year and trial factors for analysis
pred.rm.na$Year<- as.factor(pred.rm.na$Year)
pred.rm.na$Trial<-as.factor(pred.rm.na$Trial)
pred.rm.na$Treatment<-as.factor(pred.rm.na$Treatment)

# exporting data in wide format
# write.csv(pred.rm.na,"Data\\2019.1.2.predator.raw.data.csv", row.names = FALSE)

# importing egg count data for avg. number of gobies inhabiting each reef
# Will add this to predator data when I wrangle them into wide format
repro<-read.csv(file= "Data/egg.counts.2019.12.23.csv")
#making Year.fact a factor
repro$Year<-as.factor(repro$Year)
repro$Trial<-as.factor(repro$Trial)

#data wrangling####

# wide format for analyses ####

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
#View(pr)
#head(pr)

#want to drop unnecessary columns from repro df, i.e. everything but avg.inhab
pr1<-pr %>% select(Trial:Lethal.Threat,avg.inhab)
#View(pr1)

# long format for plotting #### 
# p-long = "pl"

#for all trials
pr1<-as_tibble(pr1)
pr1

pl<-pr1 %>% gather(Predator.class, Score, Present:Lethal.Threat)
#View(pl)

# Creating separate df named "HR.comp" to compare responses between HR-caged and -uncaged treatments in Trial 6
# Filtering by Trial 6 only, and removing columns for Trial and Year

HR.comp<-p %>% filter(Trial==6) %>% 
  select (-c(Trial,Year))
#View(HR.comp)

#exporting data
#write.csv(HR.comp,"Data\\2019.1.2.predator.raw.data.high.risk.comparison.csv", row.names = FALSE)

# breaking this down into two separate dataframes: one for Trials 1-5 and one for Trial 6
# Rationale: 1) only presence and sublethal predator activity is comparable in Trials 1-5
#            2) Trial 6 was testing for cage effects only, and only high-risk caged and uncaged reefs were time-lapsed

# data frame wide format to analyze Trials 1-5
t1.5.w<- pr1 %>% filter(pr1$Trial != 6)
t1.5.w$Treatment<-as.factor(t1.5.w$Treatment)
#View(t1.5.w)

# ^ this df will work for analyses for presence/absence, but not for sublethal, because no predator could ever be close enough
# to be perceived as a sublethal threat in the Low-risk treatment
# I.e. I only want to compare scores among treatments that were comparable

# Need to make a new df without Low-risk treatment in Trials 1-5
ptl.sub<-subset(t1.5.w,Treatment!="Low")
#View(ptl.sub)
ptl.sub$Treatment<-as.factor(ptl.sub$Treatment)
#levels(ptl.sub$Treatment)

#ask M. Steele: if there are differences between caged and uncaged trts
# - (e.g. presence of preds) do I have to remove that trial from the 
# - overall analyses for all trials, and run trials 1-5 and trial 6 as 
# - two separate anlayses? I think yes (see df "t1.5.l")

# will see if there are statistical differences among trials

#for trials 1-5 only (l, m, h (caged) only)
trial.1.5.long<-pl %>% filter(pl$Trial != 6)
#View(trial.1.5.long)

#for trial 6 only (HR caged and uncaged)
HR.long<-pl %>% filter(pl$Trial == 6)
#View(HR.long)

# No formal analysis will be done for lethal, but will include data in the plots to show that predators
# were seen on the high-risk reefs when they had access to them (i.e. use df "t1.5.long" for plotting)

# calculating sample sizes and number of photos that were used to calc. proportions per time lapse ####
# average number of photos that went into each proportion per reef per treatment

# a) including uncaged
#had to do it in excel pivot table...

#avg.number.of.photos
#high	11.5 (12)
#med	11.6 (12) 
#low	12.1 (12)
#uncaged	13.5 (14)

# Calculating sample sizes and arithmetic means using wide format df

# with HR caged only, using t1.4.w df, not including any data from Trial 6
Summarize(Present~Treatment,
          data=t1.5.w,
          digits=3)
# Low-risk: n = 21; mean = 0.401
# Medium-risk: n =  20; mean = 0.536
# High-risk: n =  21; mean = 0.533

# Comparing high-risk caged and uncaged treatments only in Trial 6
Summarize(Present~T6.comparison,
          data=HR.comp,
          digits=3)
# High-risk caged: n = 4; mean = 0.433 
# High-risk uncaged: n = 4; mean = 0.683

# analyses ####
#NOTE:
# analyzing mixed models with log-likelihood estimates and chi-square tests.
# start with full model, then reduce model first by non-significant random effects, then by NS. fixed effects
# want to definitely leave in Trial term as random effect, and Treatement, Year, T x Y, and avg.inhab as fixed effects.
# remove all NS. interactions with the covariate (fixed and random)

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#model structure for random effects: all slopes are the same (effects are the same among trials), but intercepts are random (magnitude differs)
## Note: use "REML=F" (maximum likelihood estimates) to compare models with log-likelihood estimates (Pinheiro & Bates, 2000; Bolker et al., 2009)

# 1 Trials 1-5, comparing predator presence and sublethal threat among Low-, Medium-, and High-risk caged treatments ####

# 1a. predator presence ####
glimpse(t1.5.w)

# full model
pre<-lmer(Present ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, t1.5.w)
hist(resid(pre))
qqnorm(resid(pre))
qqline(resid(pre))
plot(pre)
summary(pre) # 0 variance attributed to all of the random effects. 
# Will still work through log-likelihood, but if there are no effects of random variables, then I will run
# - as a regular ANCOVA, not a mixed-model ANCOVA
anova(pre)

# There may be justification for pooling in this case, or at least not including Trial in the final model
# I'm not sure if I would still do log-likelihood in that case, or if I would just run with lm

#removing three-way interaction of random effect
pre2<-update(pre, .~. -(1|avg.inhab:Treatment:Trial))
summary(pre2)
anova(pre2)

anova(pre,pre2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

pre3<-update(pre2, .~. -(1|Treatment:Trial))
summary(pre3)
anova(pre3)

anova(pre2,pre3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

pre4<-update(pre3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(pre4) # trial accounts for 6% of the residual variance
anova(pre4)

anova(pre3,pre4) # P = 0.3939, doesn't suggest any difference due to the dropping of the random effect of Trial:avg.inhab

#trial effect
pre5<-update(pre2,.~. -(1|Trial))
summary(pre5)
anova(pre5)

anova(pre2,pre5) # no sig. effect of trial

#going to pool all trials, because there is no extra variance explained by trial
#no longer dealing with a mixed model, so I have to change function to "lm"

pre6<-lm(Present ~ Treatment*Year*avg.inhab, t1.5.w)
hist(resid(pre6))
qqnorm(resid(pre6))  
qqline(resid(pre6))
summary(pre6)
anova(pre6)

anova(pre4,pre6)

# removing non-significant interactions with covariate
pre7<-update(pre6, .~. -(Treatment:Year:avg.inhab))
summary(pre7)
anova(pre7)

anova(pre7,pre6)

#now want to take out non-significant interactions with covariate (avg.inhab)

pre8<-lm(Present ~ Treatment*Year + avg.inhab, t1.5.w)
hist(resid(pre8))
qqnorm(resid(pre8))  
qqline(resid(pre8))
summary(pre8)
anova(pre8)

anova(pre7,pre8)

emmeans(pre8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(Present~Treatment,data=t1.5.w)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

# 1b. sublethal threat (compared among medium- and high-risk caged treatments) ####

# NOTE: only medium- and high-risk caged reefs in Trials 1-5 had predators that could have been perceived 
# as a sublethal threat, so I'm using the "ptl.sub" df

glimpse(ptl.sub)
levels(ptl.sub$Treatment)

# full model
st<-lmer(Sublethal.Threat ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
            (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, ptl.sub)
hist(resid(st))
qqnorm(resid(st))
qqline(resid(st))
plot(st)
summary(st) # 0 variance attributed to all of the random effects. 
# Will still work through log-likelihood, but if there are no effects of random variables, then I will run
# - as a regular ANCOVA, not a mixed-model ANCOVA
anova(st)

# There may be justification for pooling in this case, or at least not including Trial in the final model
# I'm not sure if I would still do log-likelihood in that case, or if I would just run with lm

#removing three-way interaction of random effect
st2<-update(st, .~. -(1|avg.inhab:Treatment:Trial))
summary(st2)
anova(st2)

anova(st,st2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

st3<-update(st2, .~. -(1|Treatment:Trial))
summary(st3)
anova(st3)

anova(st2,st3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

st4<-update(st3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(st4) # trial accounts for 6% of the residual variance
anova(st4)

anova(st3,st4) # P = 0.3939, doesn't suggest any difference due to the dropping of the random effect of Trial:avg.inhab

#trial effect
st5<-update(st2,.~. -(1|Trial))
summary(st5)
anova(st5)

anova(st2,st5) # no sig. effect of trial

#going to pool all trials, because there is no extra variance explained by trial
#no longer dealing with a mixed model, so I have to change function to "lm"

st6<-lm(Sublethal.Threat ~ Treatment*Year*avg.inhab, ptl.sub)
hist(resid(st6))
qqnorm(resid(st6))  
qqline(resid(st6))
summary(st6)
anova(st6)

anova(st4,st6)

# removing non-significant interactions with covariate
st7<-update(st6, .~. -(Treatment:Year:avg.inhab))
summary(st7)
anova(st7)

anova(st7,st6)

#now want to take out non-significant interactions with covariate (avg.inhab)

st8<-lm(Sublethal.Threat ~ Treatment*Year + avg.inhab, ptl.sub)
hist(resid(st8))
qqnorm(resid(st8))  
qqline(resid(st8))
summary(st8)
anova(st8)

anova(st7,st8)

emmeans(st8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(Sublethal.Threat~Treatment,data=ptl.sub)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

# 2. comparing HR caged and uncaged plots with a t-test for each distinction ####
# Rationale: Not indcluding average inhabitants because Trial 6 wsa a relatively short trial, and 
# preliminary analyses showed that the number of inhabitants did not differ than much between treatments
# (~3 gobies on average)

#loading data (in proper format for t-test) ####
hr.t.test<-read.csv(file = "Data/2020.1.30.ptl.t.test.cage.artifacts.csv")
glimpse(hr.t.test)

#renaming columns to make more sense

hr.t.test <- hr.t.test %>% 
  rename(caged_present = 1, uncaged_present = 2, caged_sublethal = 3, 
         uncaged_sublethal = 4, caged_lethal = 5, uncaged_lethal = 6)

# 2a. present
t.test(hr.t.test$caged_present,hr.t.test$uncaged_present)
# no difference in the proportion of photos with predators present between high-risk caged vs. uncaged 
# treatments

# 2b. sublethal threat
t.test(hr.t.test$caged_sublethal,hr.t.test$uncaged_sublethal)
# no diff in propotion of photos with sublethal predators

# 2c. lethal threat
t.test(hr.t.test$caged_lethal,hr.t.test$uncaged_lethal)
# no diff in lethal, in fact, equal means for each treatment

#plotting, using data in long format (trials 1-5: pl, trial6 6: HR.comp)####




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

# plotting ####

trial.1.5.long$Treatment<-ordered(trial.1.5.long$Treatment,c("Low","Medium","High"))

pred1.5<-with(trial.1.5.long, aggregate((Score),list(Treatment=Treatment,Predator.class=Predator.class),mean))
pred1.5$se<-with(trial.1.5.long, aggregate((Score),list(Treatment=Treatment,Predator.class=Predator.class), 
                                           function(x) sd(x)/sqrt(length(x))))[,3]
#ordering Predator.class values
pred1.5$Predator.class<-ordered(pred1.5$Predator.class, c("Present","Sublethal.Threat","Lethal.Threat"))
pred1.5

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
