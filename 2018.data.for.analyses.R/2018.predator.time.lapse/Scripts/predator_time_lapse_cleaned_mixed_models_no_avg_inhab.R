# Description: analyses for predator time lapses - goby reproduction MS
# Author: George C Jarvis
# Date: Tue Oct 19 14:39:15 2021
# Notes: reduced mixed models when appropriate (when model fits were singular)
# 
# Re: coding for predator activity:
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

#rm(list=ls())

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

#Column names for columns 6, 7 and 8 are not correct, need to change them to "Present", "Sublethal.Threat", 
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

# RUN ALL CODE ABOVE THIS LINE TO IMPORT AND WRANGLE DATA #######

# 1. analyses ####
#NOTE:
# analyzing mixed models with log-likelihood estimates and chi-square tests.
# start with full model, then reduce model first by non-significant random effects, then by NS. fixed effects
# want to include Trial as a random effect, and Treatement (T), Year (Y), T x Y, and avg.inhab as fixed effects.
# remove all NS. interactions with the covariate (fixed and random) at P > 0.05

# re: mixed models: when model fits were singular, did not interpret effects. Instead, removed random effect and ran as general LM
# - in general, my models were overfitted when they included all random effects and their interactions

options(contrasts = c("contr.sum","contr.poly")) #run before analyses, will set ANOVA SS to type III

#model structure for random effects: all slopes are the same (effects are the same among trials), but intercepts are random (magnitude differs)
## Note: use "REML=F" (maximum likelihood estimates) to compare models with log-likelihood estimates (Pinheiro & Bates, 2000; Bolker et al., 2009)

# 1a. Trials 1-5, comparing predator presence and sublethal threat among Low-, Medium-, and High-risk caged treatments ####

# 1ai. predator presence ####
glimpse(t1.5.w)

# pooled all trials because fit of mixed model (with trial and interactions between trial and main effects) was singular

pre6<-lm(Present ~ Treatment*Year*avg.inhab, t1.5.w)
hist(resid(pre6))
qqnorm(resid(pre6))  
qqline(resid(pre6))
#summary(pre6)
anova(pre6)

# output

# > anova(pre6)
# Analysis of Variance Table
# 
# Response: Present
#                          Df Sum Sq  Mean Sq F value Pr(>F)
# Treatment                 2 0.2469 0.123442  1.6626 0.1999
# Year                      1 0.0996 0.099582  1.3413 0.2523
# avg.inhab                 1 0.1309 0.130943  1.7637 0.1902
# Treatment:Year            2 0.0352 0.017578  0.2368 0.7901
# Treatment:avg.inhab       2 0.0388 0.019421  0.2616 0.7709
# Year:avg.inhab            1 0.0028 0.002761  0.0372 0.8479
# Treatment:Year:avg.inhab  2 0.0398 0.019897  0.2680 0.7660
# Residuals                50 3.7122 0.074244   

#anova(pre4,pre6)

# removing non-significant interactions with covariate
pre7<-update(pre6, .~. -(Treatment:Year:avg.inhab))
#summary(pre7)
anova(pre7)

# > anova(pre7)
# Analysis of Variance Table
# 
# Response: Present
# Df Sum Sq  Mean Sq F value Pr(>F)
# Treatment            2 0.2469 0.123442  1.7108 0.1907
# Year                 1 0.0996 0.099582  1.3801 0.2454
# avg.inhab            1 0.1309 0.130943  1.8148 0.1838
# Treatment:Year       2 0.0352 0.017578  0.2436 0.7847
# Treatment:avg.inhab  2 0.0388 0.019421  0.2692 0.7651
# Year:avg.inhab       1 0.0028 0.002761  0.0383 0.8457
# Residuals           52 3.7520 0.072154  

#anova(pre7,pre6)

#now want to take out non-significant effects of covariate (avg.inhab)

# (MS) final model for PTL - presence of predators, after removing NS effect of avg.inhab ####

pre8a<-lm(Present ~ Treatment*Year, t1.5.w)
# hist(resid(pre8))
# qqnorm(resid(pre8))  
# qqline(resid(pre8))
#summary(pre8)
anova(pre8a)

# Analysis of Variance Table
# 
# Response: Present
#                Df Sum Sq  Mean Sq F value Pr(>F)
# Treatment       2 0.2469 0.123442  1.7818 0.1777
# Year            1 0.0996 0.099582  1.4374 0.2356
# Treatment:Year  2 0.0802 0.040088  0.5787 0.5640
# Residuals      56 3.8795 0.069277 

# let's see if removing t*y affects the qualitative results

pre8a1<-lm(Present ~ Treatment + Year, t1.5.w)
anova(pre8a1)

# doesn't, keeping the interaction in the final model

#anova(pre7,pre8)

emmeans(pre8a, pairwise~Treatment) #warning message re: interactions, but I think it's okay

# > emmeans(pre8a, pairwise~Treatment) #warning message re: interactions, but I think it's okay
# NOTE: Results may be misleading due to involvement in interactions
# $emmeans
# Treatment emmean     SE df lower.CL upper.CL
# High       0.536 0.0636 56    0.409    0.663
# Low        0.377 0.0636 56    0.250    0.505
# Medium     0.491 0.0680 56    0.354    0.627
# 
# Results are averaged over the levels of: Year 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate     SE df t.ratio p.value
# High - Low      0.1587 0.0899 56  1.765  0.1906 
# High - Medium   0.0453 0.0931 56  0.487  0.8779 
# Low - Medium   -0.1134 0.0931 56 -1.219  0.4472 
# 
# Results are averaged over the levels of: Year 
# P value adjustment: tukey method for comparing a family of 3 estimates

boxplot(Present~Treatment,data=t1.5.w) # variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

# generating table with extracted emmeans +/-se for predator presence, based on "pre8" model

pre8a.emm <- emmeans(pre8a, ~ Treatment)

pred.pres.dfa <- as.data.frame(pre8a.emm)

# > pred.pres.dfa
# Treatment    emmean         SE df  lower.CL  upper.CL
# 1      High 0.5358690 0.06357027 56 0.4085224 0.6632155
# 2       Low 0.3771625 0.06357027 56 0.2498160 0.5045091
# 3    Medium 0.4905740 0.06795948 56 0.3544348 0.6267132

# can plot with this now, just have to do mean +/- se

pred.pres.dfa$trt_o<-ordered(pred.pres.dfa$Treatment,levels=c("Low","Medium","High"))

ptl_pres_adjusted_plot_a <- ggplot(pred.pres.dfa, aes(x=trt_o, y=emmean, fill=trt_o)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment", y = "Proportion of photos with predators present 
       (adjusted)") +
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  geom_linerange(aes(ymin=emmean-SE, ymax=emmean+SE), size=0.5,   
                 position=position_dodge(.85))

ptl_pres_adjusted_plot_a # this is the plot to plot for trials 1-5 predator presence

# START HERE NEXT TIME, NEED TO REPLOT THESE WITH AVG.INHAB REMOVED FROM MODELS FOR TRIALS 1-5 #####

# 1aii. sublethal threat (compared among medium- and high-risk caged treatments) ####

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

# MS - full model results for general linear model for sublethal threat #####

st6<-lm(Sublethal.Threat ~ Treatment*Year*avg.inhab, ptl.sub)
hist(resid(st6))
qqnorm(resid(st6))  
qqline(resid(st6))
summary(st6)
anova(st6)

# > anova(st6)
# Analysis of Variance Table
# 
# Response: Sublethal.Threat
#                          Df  Sum Sq   Mean Sq F value Pr(>F)
# Treatment                 1 0.00020 0.0002045  0.0142 0.9060
# Year                      1 0.00686 0.0068598  0.4752 0.4954
# avg.inhab                 1 0.02268 0.0226782  1.5710 0.2189
# Treatment:Year            1 0.00007 0.0000718  0.0050 0.9442
# Treatment:avg.inhab       1 0.00030 0.0002982  0.0207 0.8866
# Year:avg.inhab            1 0.00010 0.0000970  0.0067 0.9351
# Treatment:Year:avg.inhab  1 0.00009 0.0000877  0.0061 0.9384
# Residuals                33 0.47637 0.0144353

#anova(st4,st6)

# removing non-significant interactions with covariate
st7<-update(st6, .~. -(Treatment:Year:avg.inhab))
summary(st7)
anova(st7)

anova(st7,st6)

# now want to take out non-significant effects of covariate (avg.inhab)

# MS- final model for proportion of photos with sublethal predators - removed NS effect of avg.inhab ####

st8a<-lm(Sublethal.Threat ~ Treatment*Year, ptl.sub)
hist(resid(st8a))
qqnorm(resid(st8a))  
qqline(resid(st8a))
summary(st8a)
anova(st8a)

# > anova(st8a)
# Analysis of Variance Table
# 
# Response: Sublethal.Threat
#                Df  Sum Sq   Mean Sq F value Pr(>F)
# Treatment       1 0.00020 0.0002045  0.0153 0.9023
# Year            1 0.00686 0.0068598  0.5123 0.4786
# Treatment:Year  1 0.00417 0.0041713  0.3115 0.5801
# Residuals      37 0.49543 0.0133899    

emmeans(st8a, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(Sublethal.Threat~Treatment,data=ptl.sub)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

# generating table with extracted emmeans +/-se for predator presence, based on "pre8" model

st8a.emm <- emmeans(st8a, ~ Treatment)

st8a.df <- as.data.frame(st8a.emm)
st8a.df

# > st8a.df
#   Treatment    emmean         SE df   lower.CL  upper.CL
# 1      High 0.1348351 0.02794780 37 0.07820750 0.1914627
# 2    Medium 0.1186674 0.02987745 37 0.05812997 0.1792049

# can plot with this now, just have to do mean +/- se

st8a.df$trt_o<-ordered(st8a.df$Treatment,levels=c("Low","Medium","High"))

ptl_sub_adjusted_plot_a <- ggplot(st8a.df, aes(x=trt_o, y=emmean, fill=trt_o)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment", y = "Proportion of photos with predators present 
       (adjusted)") +
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  geom_linerange(aes(ymin=emmean-SE, ymax=emmean+SE), size=0.5,   
                 position=position_dodge(.85))

ptl_sub_adjusted_plot_a

# MS results - t tests for Trial 6 ######

# 1b. Trial 6, comparing HR-caged and HR-uncaged plots only with a t-test for each distinction ####
# Rationale: Not including average inhabitants because Trial 6 was a relatively short trial, and 
# preliminary analyses showed that the number of inhabitants did not differ than much between treatments
# (~3 gobies on average)

#loading data (in proper format for t-test) ####
hr.t.test<-read.csv(file = "Data/2020.1.30.ptl.t.test.cage.artifacts.csv")
glimpse(hr.t.test)

#renaming columns to make more sense

hr.t.test <- hr.t.test %>% 
  rename(caged_present = 1, uncaged_present = 2, caged_sublethal = 3, 
         uncaged_sublethal = 4, caged_lethal = 5, uncaged_lethal = 6)

# 1bi. present ####
t.test(hr.t.test$caged_present,hr.t.test$uncaged_present)
# no difference in proportion of photos with predators present between treatments: i.e. when predators were allowed full access to the reefs, 
# - (cont'd) cages did not affect the likelihood of predator presence

# > t.test(hr.t.test$caged_present,hr.t.test$uncaged_present)
# 
# Welch Two Sample t-test
# 
# data:  hr.t.test$caged_present and hr.t.test$uncaged_present
# t = -1.7617, df = 3.9878, p-value = 0.1531
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.6444856  0.1444856
# sample estimates:
#   mean of x mean of y 
# 0.4333333 0.6833333 

# 1bii. sublethal threat ####
t.test(hr.t.test$caged_sublethal,hr.t.test$uncaged_sublethal)
# no diff in propotion of photos with sublethal predators

# > t.test(hr.t.test$caged_sublethal,hr.t.test$uncaged_sublethal)
# 
# Welch Two Sample t-test
# 
# data:  hr.t.test$caged_sublethal and hr.t.test$uncaged_sublethal
# t = -1.1793, df = 4.2975, p-value = 0.2994
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.22922921  0.08994349
# sample estimates:
#   mean of x  mean of y 
# 0.04285714 0.11250000 


# 1biii. lethal threat ####
t.test(hr.t.test$caged_lethal,hr.t.test$uncaged_lethal)
# no diff in lethal, in fact, equal means for each treatment

# > t.test(hr.t.test$caged_lethal,hr.t.test$uncaged_lethal)
# 
# Welch Two Sample t-test
# 
# data:  hr.t.test$caged_lethal and hr.t.test$uncaged_lethal
# t = 4.2712e-09, df = 5.9989, p-value = 1
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1432275  0.1432275
# sample estimates:
#   mean of x  mean of y 
# 0.09166667 0.09166667 

# 2. plotting, using data in long format (trials 1-5: trial.1.5.long, trial6 6: HR.comp) ####

# 2a. Trials 1-5, including data for all treatments and predator categories ####

#ordering treatments
trial.1.5.long$Treatment<-ordered(trial.1.5.long$Treatment,c("Low","Medium","High"))

pred1.5<-with(trial.1.5.long, aggregate((Score),list(Treatment=Treatment,Predator.class=Predator.class),mean))
pred1.5$se<-with(trial.1.5.long, aggregate((Score),list(Treatment=Treatment,Predator.class=Predator.class), 
                                           function(x) sd(x)/sqrt(length(x))))[,3]
#ordering Predator.class values
pred1.5$Predator.class<-ordered(pred1.5$Predator.class, c("Present","Sublethal.Threat","Lethal.Threat"))
pred1.5

#this is the plot that I used to make the black and grayscale figure, with no figure legend, for the MS

#png("Output/2019.2.6.9.5x5.5.300dpi.png", width = 9.5, height = 5.5, units = 'in', res = 300)

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
#dev.off()

#this is the plot that I used to make a color figure, with a figure legend ("see scale fill discrete..."), for the MS
#png("Output/2019.2.6.9.5x5.5.color.300dpi.png", width = 9.5, height = 5.5, units = 'in', res = 300)

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

#dev.off()

# replotting the sublethal PTL figure for MS ####

# update: 2021/10/19 - changing values in df manually to match what they need to be in original df #####
# - want to show adjusted means (based on effects of interactions (treatment*year)  avg.inhab)
# adding means and se manually to the df, easier because plot is already made

# present
pred.pres.dfa

# > pred.pres.dfa
#   Treatment    emmean         SE df  lower.CL  upper.CL  trt_o
# 1      High 0.5358690 0.06357027 56 0.4085224 0.6632155   High
# 2       Low 0.3771625 0.06357027 56 0.2498160 0.5045091    Low
# 3    Medium 0.4905740 0.06795948 56 0.3544348 0.6267132 Medium

# sublethal threat
st8a.df

# > st8a.df
#   Treatment    emmean         SE df   lower.CL  upper.CL  trt_o
# 1      High 0.1348351 0.02794780 37 0.07820750 0.1914627   High
# 2    Medium 0.1186674 0.02987745 37 0.05812997 0.1792049 Medium

# changing values in 'pred.6' df (it's formulated below) 

# 1. mean
# 2. se

# present

# low 
pred1.5[4, 3] = 0.377
pred1.5[4, 4] = 0.06
# med
pred1.5[5, 3] = 0.491
pred1.5[5, 4] = 0.07
# high
pred1.5[6, 3] = 0.536
pred1.5[6, 4] = 0.06

# sublethal

# medium
pred1.5[8, 3] = 0.117
pred1.5[8, 4] = 0.03
# high
pred1.5[9, 3] = 0.135
pred1.5[9, 4] = 0.03

# replotting

# for grayscale plot (no legend)

png("Output/PTL_trials_1_5_LS_means_no_avg_inhab.png", width = 9.5, height = 5.5, units = 'in', res = 300)

t1.5.plot<- ggplot(pred1.5, aes(x=Treatment, y=x, fill=Predator.class)) +
  geom_bar(stat="identity", colour= "black", width = 0.9, position="dodge", show.legend = FALSE)+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment", y = "Proportion of photos with predators present 
       (adjusted)")+ 
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

# to get legend as you want it (but don't need to repeat this, because already done)

#png("Output/2019.2.6.9.5x5.5.color.300dpi.png", width = 9.5, height = 5.5, units = 'in', res = 300)

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

#dev.off()

# 2b. Trial 6 plots comparing high-risk caged and uncaged treatments (HR.long df) ####
# NOTE: Using "T6.comparison" for treatment factor this time

#View(HR.long)

HR.long$T6.comparison<-as.factor(HR.long$T6.comparison)

pred.6<-with(HR.long, aggregate((Score),list(T6.comparison=T6.comparison,Predator.class=Predator.class),mean))
pred.6$se<-with(HR.long, aggregate((Score),list(T6.comparison=T6.comparison,Predator.class=Predator.class), 
                                   function(x) sd(x)/sqrt(length(x))))[,3]
pred.6

pred.6$Predator.class<-as.factor(pred.6$Predator.class)
#ordering Predator.class values
pred.6$Predator.class<-ordered(pred.6$Predator.class, c("Present","Sublethal.Threat","Lethal.Threat"))

pred.6<- fct_relevel(pred.6$Predator.class,"Present","Sublethal.Threat","Lethal.Threat")
pred.6

#renaming the labels for high-risk treatment t6.comparison with base R from High to "High - Caged" and "Uncaged" to "High - Uncaged" 
#will eventually have a second figure that has "Risk Treatment" as the x axis, and can just compare caging effects, that will be more clear anyway
#doing it in the new df ("pred.6"), not in the original df used to make calculations ("HR.long")

# Rename by name: change "High" to "Caged"
levels(pred.6$T6.comparison)[levels(pred.6$T6.comparison)=="High"] <- "High - Caged"
levels(pred.6$T6.comparison)[levels(pred.6$T6.comparison)=="Uncaged"] <- "High - Uncaged"

#png("Output/2019.2.6.9.5x5.5.trial6.300dpi.png", width = 9.5, height = 5.5, units = 'in', res = 300)

t6.plot<- ggplot(pred.6, aes(x=T6.comparison, y=x, fill=Predator.class)) +
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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.85))
t6.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.60)) + theme(text = element_text(family="Arial"))
#dev.off()
