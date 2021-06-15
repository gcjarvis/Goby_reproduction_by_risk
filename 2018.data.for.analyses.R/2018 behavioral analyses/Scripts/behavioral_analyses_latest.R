# Description: Script for behavioral analyses from Jarvis and Steele 
# Author: George C Jarvis
# Date: Tue Jun 15 17:03:24 2021
# Notes: No behavioral observations for Trial 6, so no comparisons of behaviors between high-risk caged and uncaged treatments
# --------------

# Description: Script for behavioral analyses from Jarvis and Steele 
# Author: George C Jarvis
# Date: Sat Apr 11 17:25:38 2020
# Notes: No behavioral observations for Trial 6, so no comparisons of behaviors between high-risk
## caged and uncaged treatments
# --------------

# load packages

library(lme4)
# library(lmerTest) # I think I should go with analysis of deviance for mixed model ANCOVA's
library(car)
library(tidyverse)
library(emmeans) #for generating least-squares adjusted means from models 
library(ggpubr) # version 0.4.0
library(extrafont)
#extrafont::font_import()

#importing data ####
behave<-read.csv("Data/2019.10.25.behavior.includes.recollections.csv")

# data manipulation ####
# adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
# NOTE: there were no behavioral observations for high-risk uncaged reefs in trial 6
behave$Year <- ifelse(behave$Trial <=3, 2017, 2018)
#making Year and Trial factors
behave$Year<- as.factor(behave$Year)
behave$Trial<- as.factor(behave$Trial)

#making the variable "avg.inhab" ((20+reco)/2), rounded to the nearest whole fish
#using the average number of inhabitants per reef as the covariate in mixed models
behave$avg.inhab<-(ceiling((behave$Recollection+20)/2))

# 1. analyses ####

# analyzing with mixed mixed model ANCOVA's, with Trial as a random effect 

# tested random effect with with likelihood ratio tests, and fixed effects were tested with analysis
# start with full model, then reduce model first by non-significant random effects, then by NS. fixed effects
# At the very least, I will include Trial term as random effect, and Treatment, Year, T x Y, and avg.inhab as fixed effects.
# remove all NS. interactions with the covariate (fixed and random)

#NOTE: the values for estimates that I included in the results (supplementary table for behaviors)
## - came from the summary estimates in the fully-reduced model for fixed factors; 
## - (random factors don't have estimates)

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

# 1a. proportion of time exposed ####

# full model with all possible interactions between fixed and random effects (fit is singular)
pe<-lmer(proportion.exposed~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(pe))
qqnorm(resid(pe))
qqline(resid(pe))
plot(pe)
summary(pe) #don't seem to be very many differences based on trial
#anova(pe)
Anova(pe, type = "III")

rand(pe) # none of the residual variance in exposure time is explained by the random effects, running as fixed effects ANCOVA

# model selection, minus the interactions between fixed and random effects#

pe.1<-lmer(proportion.exposed~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab) + (Treatment*Year*avg.inhab)+
           (1|Trial), data = behave)
summary(pe.1)
Anova(pe.1, type = "III")

# result:

# > Anova(pe.1, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: proportion.exposed
# Chisq Df Pr(>Chisq)   
# (Intercept)              7.9097  1   0.004917 **
#   Treatment                0.0781  2   0.961715   
# Year                     0.3830  1   0.535989   
# avg.inhab                0.8308  1   0.362050   
# Treatment:Year           1.9122  2   0.384385   
# Treatment:avg.inhab      0.3940  2   0.821194   
# Year:avg.inhab           0.2493  1   0.617561   
# Treatment:Year:avg.inhab 1.8093  2   0.404682   

# remove 3-way 

pe.1a<-lmer(proportion.exposed~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab) +
             (1|Trial), data = behave)
summary(pe.1a)
Anova(pe.1a, type = "III")

# result

# > Anova(pe.1a, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: proportion.exposed
# Chisq Df Pr(>Chisq)   
# (Intercept)         8.3433  1   0.003871 **
#   Treatment           0.1363  2   0.934132   
# Year                0.5199  1   0.470896   
# avg.inhab           0.7770  1   0.378068   
# Treatment:Year      0.4946  2   0.780911   
# Treatment:avg.inhab 0.1110  2   0.946010   
# Year:avg.inhab      0.3444  1   0.557316 

# remove year*avg.inhab 

pe.1b<-lmer(proportion.exposed~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) +
              (1|Trial), data = behave)
summary(pe.1b)
Anova(pe.1b, type = "III")

# result:

# > Anova(pe.1b, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: proportion.exposed
# Chisq Df Pr(>Chisq)    
# (Intercept)         12.0203  1  0.0005262 ***
#   Treatment            0.1705  2  0.9183021    
# Year                 1.0447  1  0.3067312    
# avg.inhab            0.3760  1  0.5397304    
# Treatment:Year       0.5688  2  0.7524569    
# Treatment:avg.inhab  0.1574  2  0.9243358   

# remove treatment*avg.inhab 

pe.1c<-lmer(proportion.exposed~ Treatment + Year + avg.inhab + (Treatment*Year) +
              (1|Trial), data = behave)
summary(pe.1c)
Anova(pe.1c, type = "III")

# result:

# > Anova(pe.1c, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: proportion.exposed
# Chisq Df Pr(>Chisq)    
# (Intercept)    12.3602  1  0.0004386 ***
# Treatment       8.6600  2  0.0131675 *  
# Year            1.0629  1  0.3025598    
# avg.inhab       0.5897  1  0.4425333    
# Treatment:Year  0.4941  2  0.7811050  

# remove treatment*year (this is final model)

pe.1d<-lmer(proportion.exposed~ Treatment + Year + avg.inhab +
              (1|Trial),  data = behave)
summary(pe.1d) # trial explain 3% of the residual variance, low, but think it's worth including
Anova(pe.1d, type = "III")

# result:

# > Anova(pe.1d, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: proportion.exposed
# Chisq Df Pr(>Chisq)    
# (Intercept) 13.2363  1  0.0002746 ***
#   Treatment    9.1129  2  0.0104992 *  
#   Year         1.0804  1  0.2986023    
#   avg.inhab    0.5180  1  0.4716895    

# go with these results for in-text reporting

# testing effect of trial with likelihood ratio test

# making dummy variable with randomized values for trial

set.seed(595)
behave$trial_rand <- sample(behave$Trial)

# running model with trial and dummy variable included as random effects

pe.1d_dr<-lmer(proportion.exposed~ Treatment + Year + avg.inhab +
              (1|Trial) + (1|trial_rand), data = behave)
summary(pe.1d_dr) # trial explain 3% of the residual variance, low, but think it's worth including
Anova(pe.1d_dr, type = "III")

# removing trial effect, then running likelihood ratio test

pe.1d_nt<-lmer(proportion.exposed~ Treatment + Year + avg.inhab + (1|trial_rand), data = behave)
summary(pe.1d_nt)

2*(logLik(pe.1d_dr) - logLik(pe.1d_nt)) # Chi2 =  0.27
pchisq(2*(logLik(pe.1d_dr) - logLik(pe.1d_nt)), df = 1, lower.tail=F) # P = 0.60 df = 1

# trial not significant, but this looks like final model for now

# looks like final model, plotting this with emmeans to see how LS-means for treatment compare to plots of raw data

emmeans(pe.1d, pairwise ~ Treatment)

# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# > emmeans(pe.1d, pairwise ~ Treatment)
# $emmeans
# Treatment emmean     SE   df lower.CL upper.CL
# High       0.585 0.0388 16.1    0.502    0.667
# Low        0.731 0.0388 16.1    0.649    0.813
# Medium     0.687 0.0386 17.4    0.606    0.768
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate     SE   df t.ratio p.value
# High - Low     -0.1462 0.0497 84.3 -2.941  0.0116 
# High - Medium  -0.1023 0.0499 84.8 -2.051  0.1064 
# Low - Medium    0.0439 0.0499 84.9  0.879  0.6548 
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# plotting

pe.1d.emm <- emmeans(pe.1d, ~ Treatment)
plot(pe.1d.emm)

# LS - means vs. raw data for proportion of time exposed

bargraph.CI(x.factor = Treatment, response = proportion.exposed, main="raw data - proportion exposed vs. treatment", 
            xlab="Treatment", ylab="proportion of time exposed", data = behave)

# raw data (i.e. unadjusted plots) are fine, plot raw data for exposure behavior

# 1b. linear distance traveled ####

dm<-lmer(total.dist.moved~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(dm))
qqnorm(resid(dm))
qqline(resid(dm))
plot(dm)
boxplot(total.dist.moved~Treatment,data=behave) #slightly higher variance in Low treatment
summary(dm) #seems to be a bit more variation by trial
anova(dm)

rand(dm) #trial is significant, all others are not, so going to remove them from model

# model selection

dm.1a<-lmer(total.dist.moved~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab) + (Treatment*Year*avg.inhab) +
              (1|Trial), REML=F, data = behave)
Anova(dm.1a, type = "III")

# result:

# > Anova(dm.1a, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: total.dist.moved
# Chisq Df Pr(>Chisq)  
# (Intercept)              0.0593  1    0.80763  
# Treatment                0.2668  2    0.87510  
# Year                     0.3525  1    0.55270  
# avg.inhab                6.3378  1    0.01182 *
#   Treatment:Year           6.6488  2    0.03599 *
#   Treatment:avg.inhab      0.8185  2    0.66416  
# Year:avg.inhab           0.5030  1    0.47818  
# Treatment:Year:avg.inhab 5.5622  2    0.06197 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# removing 3-way

dm.1b<-lmer(total.dist.moved~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab) +
              (1|Trial), REML=F, data = behave)
Anova(dm.1b, type = "III")

# result:

# > Anova(dm.1b, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: total.dist.moved
# Chisq Df Pr(>Chisq)  
# (Intercept)         0.0019  1     0.9649  
# Treatment           0.1725  2     0.9173  
# Year                0.8785  1     0.3486  
# avg.inhab           4.7677  1     0.0290 *
#   Treatment:Year      2.7677  2     0.2506  
# Treatment:avg.inhab 0.0474  2     0.9766  
# Year:avg.inhab      1.2222  1     0.2689  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# removing year*avg.inhab

dm.1c<-lmer(total.dist.moved~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) +
              (1|Trial), REML=F, data = behave)
Anova(dm.1c, type = "III")

# result:

# > Anova(dm.1c, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: total.dist.moved
# Chisq Df Pr(>Chisq)   
# (Intercept)         0.1781  1   0.672969   
# Treatment           0.1545  2   0.925673   
# Year                0.1190  1   0.730070   
# avg.inhab           8.5815  1   0.003396 **
#   Treatment:Year      2.3974  2   0.301582   
# Treatment:avg.inhab 0.0852  2   0.958277   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# removing trt*avg.inhab

dm.1d<-lmer(total.dist.moved~ Treatment + Year + avg.inhab + (Treatment*Year) +
              (1|Trial), REML=F, data = behave)
Anova(dm.1d, type = "III")

# result:

# > Anova(dm.1d, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: total.dist.moved
# Chisq Df Pr(>Chisq)   
# (Intercept)    0.1318  1    0.71658   
# Treatment      6.3134  2    0.04257 * 
#   Year           0.1170  1    0.73233   
# avg.inhab      8.7835  1    0.00304 **
#   Treatment:Year 2.7251  2    0.25601   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# removing trt*year

dm.1e<-lmer(total.dist.moved~ Treatment + Year + avg.inhab +
              (1|Trial), REML=F, data = behave)
Anova(dm.1e, type = "III")

# result:

# > Anova(dm.1e, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: total.dist.moved
# Chisq Df Pr(>Chisq)   
# (Intercept) 0.0756  1   0.783309   
# Treatment   5.1573  2   0.075876 . 
# Year        0.1004  1   0.751298   
# avg.inhab   8.1530  1   0.004299 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# trying same model, but with REML = T

dm.1e_remlt<-lmer(total.dist.moved~ Treatment + Year + avg.inhab +
              (1|Trial), REML=T, data = behave)
Anova(dm.1e_remlt, type = "III")

# result:

# > Anova(dm.1e_remlt, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: total.dist.moved
# Chisq Df Pr(>Chisq)   
# (Intercept) 0.0662  1   0.796961   
# Treatment   4.9348  2   0.084805 . 
# Year        0.0331  1   0.855729   
# avg.inhab   7.6644  1   0.005632 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# same answer, go with this model

# estimate for avg inhab?

summary(dm.1e_remlt)

# result:

# > summary(dm.1e_remlt)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: total.dist.moved ~ Treatment + Year + avg.inhab + (1 | Trial)
#    Data: behave
# 
# REML criterion at convergence: 1030.5
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -1.43379 -0.84374 -0.07003  0.70481  2.54237 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Trial    (Intercept) 3900     62.45   
#  Residual             4936     70.26   
# Number of obs: 93, groups:  Trial, 6
# 
# Fixed effects:
#             Estimate Std. Error      df t value Pr(>|t|)   
# (Intercept)  -16.039     62.342  55.182  -0.257  0.79792   
# Treatment1   -11.648     10.315  83.388  -1.129  0.26206   
# Treatment2    22.941     10.328  83.404   2.221  0.02905 * 
# Year1         -4.881     26.848   3.373  -0.182  0.86612   
# avg.inhab     12.075      4.361  85.305   2.768  0.00691 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#            (Intr) Trtmn1 Trtmn2 Year1 
# Treatment1 -0.028                     
# Treatment2 -0.052 -0.495              
# Year1      -0.054 -0.004 -0.002       
# avg.inhab  -0.903  0.033  0.060  0.044

# linear distance traveled was greater when there were more gobies on the reef - neat!

# looking at LS-means for estimates of treatment

emmeans(dm.1e_remlt, pairwise~Treatment)

# results:

# > emmeans(dm.1e_remlt, pairwise~Treatment)
# $emmeans
# Treatment emmean   SE   df lower.CL upper.CL
# High         129 28.8 5.22     55.3      202
# Low          163 28.9 5.23     89.9      236
# Medium       129 28.7 5.16     55.8      202
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate   SE   df t.ratio p.value
# High - Low     -34.589 17.8 84.0 -1.938  0.1344 
# High - Medium   -0.355 17.9 84.2 -0.020  0.9998 
# Low - Medium    34.235 18.0 84.2  1.906  0.1433 
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# plot:

dm.1e_remlt.emm <- emmeans(dm.1e_remlt, ~ Treatment)
plot(dm.1e_remlt.emm) # slightly different from raw means, and sig. effect of avg.inhab, so I think we go with adjusted means

# making df

dm.1e.emm <- emmeans(dm.1e_remlt, pairwise ~ Treatment) # to be stored as df
dm.1e.emm.df <- as.data.frame(dm.1e.emm$emmeans)
View(dm.1e.emm.df)

# testing effect of trial with likelihood ratio test

# making dummy variable with randomized values for trial

set.seed(595)
behave$trial_rand <- sample(behave$Trial)

# running model with trial and dummy variable included as random effects

dm.1e_remlt_dr<-lmer(total.dist.moved~ Treatment + Year + avg.inhab +
                 (1|Trial) + (1|trial_rand), data = behave)
summary(dm.1e_remlt_dr) # trial explain 3% of the residual variance, low, but think it's worth including
Anova(dm.1e_remlt_dr, type = "III")

# removing trial effect, then running likelihood ratio test

dm.1e_remlt_nt<-lmer(total.dist.moved~ Treatment + Year + avg.inhab + (1|trial_rand), data = behave)
summary(dm.1e_remlt_nt)

2*(logLik(dm.1e_remlt_dr) - logLik(dm.1e_remlt_nt)) # Chi2 =  13.69
pchisq(2*(logLik(dm.1e_remlt_dr) - logLik(dm.1e_remlt_nt)), df = 1, lower.tail=F) # P < 0.01 df = 1

ranef(dm.1e_remlt) 

# results: this, plus likelihood ratio tests for random effect of trial suggest differences in diatsnce moved among trials

# > ranef(dm.1e_remlt)
# $Trial
# (Intercept)
# 1    4.974796
# 2  -11.201416
# 3    6.226620
# 4    4.332533
# 5   80.964289
# 6  -85.296822
# 
# with conditional variances for “Trial”

# 1c. foraging rate ####
fr<-lmer(bites.min~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(fr))
qqnorm(resid(fr))
qqline(resid(fr))
plot(fr)
boxplot(bites.min~Treatment,data=behave) #slightly higher variance in Low treatment
summary(fr) #seems to be a bit more variation by trial
anova(fr)

rand(fr) #trial explains some of the variance, but not much, fit turns out to be singluar if I try to run as mixed ANCOVA

# model selection as ficed ANCOVA instead:

# full model

fr_fix <-lm(bites.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab) + (Treatment*Year*avg.inhab), data = behave)
Anova(fr_fix, type = "III")

# results:

# > Anova(fr_fix, type = "III")
# Anova Table (Type III tests)
# 
# Response: bites.min
#                           Sum Sq Df F value Pr(>F)
# (Intercept)               0.4280  1  1.5106 0.2226
# Treatment                 0.5555  2  0.9804 0.3796
# Year                      0.3803  1  1.3423 0.2500
# avg.inhab                 0.0287  1  0.1013 0.7511
# Treatment:Year            0.0508  2  0.0897 0.9143
# Treatment:avg.inhab       0.4700  2  0.8295 0.4400
# Year:avg.inhab            0.3867  1  1.3649 0.2461
# Treatment:Year:avg.inhab  0.0174  2  0.0307 0.9697
# Residuals                22.9497 81   

# rem three-way

fr_fix1 <-lm(bites.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab), data = behave)
Anova(fr_fix1, type = "III")

# result:

# > Anova(fr_fix1, type = "III")
# Anova Table (Type III tests)
# 
# Response: bites.min
# Sum Sq Df F value Pr(>F)
# (Intercept)          0.4424  1  1.5986 0.2096
# Treatment            0.5785  2  1.0453 0.3562
# Year                 0.4391  1  1.5870 0.2113
# avg.inhab            0.0291  1  0.1050 0.7467
# Treatment:Year       0.7169  2  1.2955 0.2793
# Treatment:avg.inhab  0.5018  2  0.9068 0.4078
# Year:avg.inhab       0.4479  1  1.6188 0.2068
# Residuals           22.9672 83 

# rem year*avg.inhab

fr_fix2 <-lm(bites.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab), data = behave)
Anova(fr_fix2, type = "III")

# result:

# > Anova(fr_fix2, type = "III")
# Anova Table (Type III tests)
# 
# Response: bites.min
# Sum Sq Df F value  Pr(>F)  
# (Intercept)          0.8277  1  2.9694 0.08853 .
# Treatment            0.6845  2  1.2277 0.29816  
# Year                 0.0000  1  0.0001 0.99167  
# avg.inhab            0.0028  1  0.0100 0.92061  
# Treatment:Year       0.8375  2  1.5023 0.22854  
# Treatment:avg.inhab  0.5943  2  1.0660 0.34901  
# Residuals           23.4151 84                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# rem trt*avg.inhab

fr_fix3 <-lm(bites.min ~ Treatment + Year + avg.inhab + (Treatment*Year), data = behave)
Anova(fr_fix3, type = "III")

# result:

# > Anova(fr_fix3, type = "III")
# Anova Table (Type III tests)
# 
# Response: bites.min
# Sum Sq Df F value Pr(>F)
# (Intercept)     0.5643  1  2.0213 0.1587
# Treatment       0.6094  2  1.0914 0.3404
# Year            0.0021  1  0.0075 0.9312
# avg.inhab       0.0225  1  0.0804 0.7774
# Treatment:Year  0.6355  2  1.1382 0.3252
# Residuals      24.0094 86 

# rem trt*year

fr_fix4 <-lm(bites.min ~ Treatment + Year + avg.inhab, data = behave)
Anova(fr_fix4, type = "III")

# results (final model):

# > Anova(fr_fix4, type = "III")
# Anova Table (Type III tests)
# 
# Response: bites.min
# Sum Sq Df F value Pr(>F)
# (Intercept)  0.5939  1  2.1207 0.1489
# Treatment    0.6946  2  1.2402 0.2943
# Year         0.0023  1  0.0084 0.9274
# avg.inhab    0.0191  1  0.0682 0.7946
# Residuals   24.6449 88   

# emmeans

emmeans(fr_fix4, pairwise~Treatment)

# results:

# > emmeans(fr_fix4, pairwise~Treatment)
# $emmeans
# Treatment emmean     SE df lower.CL upper.CL
# High       0.814 0.0955 88    0.624    1.004
# Low        0.613 0.0955 88    0.423    0.803
# Medium     0.656 0.0958 88    0.465    0.846
# 
# Results are averaged over the levels of: Year 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate    SE df t.ratio p.value
# High - Low      0.2010 0.134 88  1.495  0.2982 
# High - Medium   0.1581 0.135 88  1.174  0.4717 
# Low - Medium   -0.0429 0.135 88 -0.318  0.9458 
# 
# Results are averaged over the levels of: Year 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# plot this next and see how it compares to plots fo the raw data

fr_fix4.emm <- emmeans(fr_fix4, ~ Treatment)
plot(fr_fix4.emm)

# storing as df, I think plots of raw data might be misleading

fr_fix4.emm.pw <- emmeans(fr_fix4, pairwise ~ Treatment) # to be stored as df
fr_fix4.emm.df <- as.data.frame(fr_fix4.emm.pw$emmeans)
View(fr_fix4.emm.df)

# as a ggplot (I think the raw data are just fine)

fr<- ggplot(fr_fix4.emm.df, aes(x=Treatment, y=emmean, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment",y=(expression(atop("Foraging rate"~(bites~min^-1), 
                                             paste("(adjusted)")))))+
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
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.01)) +
  geom_linerange(aes(ymin=emmean-SE, ymax=emmean+SE), size=0.5,   
                 position=position_dodge(.85))# + theme(text = element_text(family="Arial"))


fr

# 1d. interactins with conspecifics ####

ci<-lmer(courtship.min~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(ci))
qqnorm(resid(ci))
qqline(resid(ci))
plot(ci)
boxplot(courtship.min~Treatment,data=behave) #slightly higher variance in Low treatment
summary(ci) #seems to be a bit more variation by trial
anova(ci)

rand(ci) #will keep trial, explains some of the variance

# model selection, including trial as random effect

# full model

ci.a <-lmer(courtship.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab) + (Treatment*Year*avg.inhab) +
              (1|Trial), REML=F, data = behave)
Anova(ci.a, type = "III")

# results:

# > Anova(ci.a, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: courtship.min
# Chisq Df Pr(>Chisq)  
# (Intercept)              4.2804  1    0.03855 *
#   Treatment                3.5651  2    0.16821  
# Year                     1.4450  1    0.22934  
# avg.inhab                0.4731  1    0.49157  
# Treatment:Year           3.6158  2    0.16400  
# Treatment:avg.inhab      3.4157  2    0.18126  
# Year:avg.inhab           1.2097  1    0.27139  
# Treatment:Year:avg.inhab 3.4089  2    0.18187  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# remove 3-way

ci.a1 <-lmer(courtship.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab) +
              (1|Trial), REML=F, data = behave)
Anova(ci.a1, type = "III")

# results:

# > Anova(ci.a1, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: courtship.min
# Chisq Df Pr(>Chisq)  
# (Intercept)         5.6811  1    0.01715 *
#   Treatment           2.8592  2    0.23940  
# Year                0.5725  1    0.44927  
# avg.inhab           0.9923  1    0.31918  
# Treatment:Year      0.6994  2    0.70490  
# Treatment:avg.inhab 2.8464  2    0.24094  
# Year:avg.inhab      0.3954  1    0.52945  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# remove year*avg.inhab

ci.a2 <-lmer(courtship.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) +
               (1|Trial), REML=F, data = behave)
Anova(ci.a2, type = "III")

# results:

# > Anova(ci.a2, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: courtship.min
# Chisq Df Pr(>Chisq)   
# (Intercept)         7.4785  1   0.006244 **
#   Treatment           2.7298  2   0.255410   
# Year                0.7648  1   0.381817   
# avg.inhab           1.7275  1   0.188728   
# Treatment:Year      0.7206  2   0.697473   
# Treatment:avg.inhab 2.7154  2   0.257251   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# remove trt*avg.inhab

ci.a3 <-lmer(courtship.min ~ Treatment + Year + avg.inhab + (Treatment*Year) +
               (1|Trial), REML=F, data = behave)
Anova(ci.a3, type = "III")

# results:

# > Anova(ci.a3, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: courtship.min
# Chisq Df Pr(>Chisq)   
# (Intercept)    6.7142  1   0.009565 **
#   Treatment      3.4617  2   0.177132   
# Year           0.8798  1   0.348255   
# avg.inhab      1.3070  1   0.252932   
# Treatment:Year 0.7587  2   0.684292   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# remove trt*year

ci.a4 <-lmer(courtship.min ~ Treatment + Year + avg.inhab +
               (1|Trial), REML=T, data = behave)
Anova(ci.a4, type = "III")

# results:

# > Anova(ci.a4, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: courtship.min
# Chisq Df Pr(>Chisq)   
# (Intercept) 7.1933  1   0.007318 **
#   Treatment   4.0457  2   0.132276   
# Year        0.4498  1   0.502443   
# avg.inhab   1.7158  1   0.190238   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# testing effect of trial with likelihood ratio test

# making dummy variable with randomized values for trial

set.seed(595)
behave$trial_rand <- sample(behave$Trial)

# running model with trial and dummy variable included as random effects

ci.a4_dr<-lmer(courtship.min~ Treatment + Year + avg.inhab +
                       (1|Trial) + (1|trial_rand), data = behave)
summary(ci.a4_dr) # trial explain 3% of the residual variance, low, but think it's worth including
Anova(ci.a4_dr, type = "III")

# removing trial effect, then running likelihood ratio test

ci.a4_dr_nt<-lmer(courtship.min~ Treatment + Year + avg.inhab + (1|trial_rand), data = behave)
summary(ci.a4_dr_nt)

2*(logLik(ci.a4_dr) - logLik(ci.a4_dr_nt)) # Chi2 =  1.18
pchisq(2*(logLik(ci.a4_dr) - logLik(ci.a4_dr_nt)), df = 1, lower.tail=F) # P =0.28 df = 1

ranef(ci.a4) 

# results (all fairly similar among trials):

# > ranef(ci.a4) 
# $Trial
# (Intercept)
# 1  0.011421816
# 2 -0.020170061
# 3  0.008748245
# 4 -0.004974428
# 5  0.015663868
# 6 -0.010689440
# 
# with conditional variances for “Trial”

# looking at LS-means

emmeans(ci.a4, pairwise~Treatment)

# results:

# > emmeans(ci.a4, pairwise~Treatment)
# $emmeans
# Treatment emmean     SE   df lower.CL upper.CL
# High      0.0799 0.0189 45.2   0.0419    0.118
# Low       0.0939 0.0189 45.2   0.0559    0.132
# Medium    0.1237 0.0187 47.0   0.0860    0.161
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate     SE   df t.ratio p.value
# High - Low     -0.0140 0.0223 89.9 -0.626  0.8063 
# High - Medium  -0.0437 0.0224 90.6 -1.955  0.1295 
# Low - Medium   -0.0298 0.0224 90.7 -1.330  0.3825 
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# plot this next and see how it compares to plots fo the raw data

ci.a4.emm <- emmeans(ci.a4, ~ Treatment)
plot(ci.a4.emm)

# looks very similar to the plots of raw data for courtship rates, go with raw data

# 1e. movement rate ####
mr<-lmer(movements.min~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(mr))
qqnorm(resid(mr))
qqline(resid(mr))
plot(mr)
boxplot(movements.min~Treatment,data=behave) #slightly higher variance in Low treatment
summary(mr) #seems to be a bit more variation by trial
anova(mr)

rand(mr) #going to remove trial, seems like it doesn't explain much variation
# P = 0.94

# mixed ANCOVA has singular fit, running as fixed ANCOVA instead

mr.a_f <-lm(movements.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab) + (Treatment*Year*avg.inhab), data = behave)
Anova(mr.a_f, type = "III")

# results:

# > Anova(mr.a_f, type = "III")
# Anova Table (Type III tests)
# 
# Response: movements.min
# Sum Sq Df F value  Pr(>F)  
# (Intercept)               0.9607  1  3.1480 0.07978 .
# Treatment                 0.0094  2  0.0155 0.98466  
# Year                      0.3422  1  1.1215 0.29274  
# avg.inhab                 0.1782  1  0.5840 0.44697  
# Treatment:Year            0.0998  2  0.1636 0.84938  
# Treatment:avg.inhab       0.0034  2  0.0056 0.99445  
# Year:avg.inhab            0.3014  1  0.9875 0.32331  
# Treatment:Year:avg.inhab  0.0881  2  0.1443 0.86583  
# Residuals                24.7190 81                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# remove 3-way

mr.a_f1 <-lm(movements.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) + (Year*avg.inhab), data = behave)
Anova(mr.a_f1, type = "III")

# results:

# > Anova(mr.a_f1, type = "III")
# Anova Table (Type III tests)
# 
# Response: movements.min
# Sum Sq Df F value  Pr(>F)  
# (Intercept)          1.0950  1  3.6635 0.05906 .
# Treatment            0.0377  2  0.0631 0.93885  
# Year                 0.2964  1  0.9916 0.32224  
# avg.inhab            0.1431  1  0.4789 0.49086  
# Treatment:Year       0.1108  2  0.1854 0.83110  
# Treatment:avg.inhab  0.0072  2  0.0120 0.98806  
# Year:avg.inhab       0.2541  1  0.8500 0.35922  
# Residuals           24.8071 83     

# remove year*avg.inhab

mr.a_f2 <-lm(movements.min ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab), data = behave)
#summary(mr.a_f2)
Anova(mr.a_f2, type = "III")

# results:

# Anova Table (Type III tests)
# 
# Response: movements.min
# Sum Sq Df F value  Pr(>F)  
# (Intercept)          1.5811  1  5.2994 0.02381 *
#   Treatment            0.0640  2  0.1072 0.89843  
# Year                 0.1148  1  0.3848 0.53673  
# avg.inhab            0.0506  1  0.1696 0.68153  
# Treatment:Year       0.1224  2  0.2051 0.81501  
# Treatment:avg.inhab  0.0185  2  0.0311 0.96943  
# Residuals           25.0611 84

# remove trt*avg.inhab

mr.a_f3 <-lm(movements.min ~ Treatment + Year + avg.inhab + (Treatment*Year), data = behave)
Anova(mr.a_f3, type = "III")

# results:

# > Anova(mr.a_f3, type = "III")
# Anova Table (Type III tests)
# 
# Response: movements.min
# Sum Sq Df F value  Pr(>F)  
# (Intercept)     1.6321  1  5.5964 0.02025 *
#   Treatment       1.1307  2  1.9385 0.15015  
# Year            0.1203  1  0.4124 0.52244  
# avg.inhab       0.0701  1  0.2404 0.62516  
# Treatment:Year  0.1130  2  0.1938 0.82421  
# Residuals      25.0797 86 

# remove trt*year (final model)

mr.a_f4 <-lm(movements.min ~ Treatment + Year + avg.inhab, data = behave)
Anova(mr.a_f4, type = "III")

# results:

# > Anova(mr.a_f4, type = "III")
# Anova Table (Type III tests)
# 
# Response: movements.min
# Sum Sq Df F value  Pr(>F)  
# (Intercept)  1.7771  1  6.2077 0.01459 *
#   Treatment    1.1883  2  2.0754 0.13161  
# Year         0.1266  1  0.4423 0.50775  
# avg.inhab    0.0494  1  0.1726 0.67884  
# Residuals   25.1927 88  

# looking at LS-means

emmeans(mr.a_f4, pairwise~Treatment)

# results:

# $emmeans
# Treatment emmean     SE df lower.CL upper.CL
# High        1.10 0.0965 88    0.910     1.29
# Low         1.35 0.0965 88    1.156     1.54
# Medium      1.11 0.0968 88    0.922     1.31
# 
# Results are averaged over the levels of: Year 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate    SE df t.ratio p.value
# High - Low     -0.2463 0.136 88 -1.812  0.1717 
# High - Medium  -0.0129 0.136 88 -0.095  0.9951 
# Low - Medium    0.2334 0.136 88  1.712  0.2063 
# 
# Results are averaged over the levels of: Year 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# plot this next and see how it compares to plots fo the raw data

mr.a_f4.emm <- emmeans(mr.a_f4, ~ Treatment)
plot(mr.a_f4.emm)

# looks very similar to the plots of raw data for courtship rates, plot raw data

# 2. plotting ####

# proportion of time exposed ####
exp<-with(behave, aggregate((proportion.exposed), list(Treatment=Treatment), mean))
exp$se<-with(behave, aggregate((proportion.exposed), list(Treatment=Treatment), 
                               function(x) sd(x)/sqrt(length(x))))[,2]

e<- ggplot(exp, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment", y="Proportion of time exposed") +
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
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.807)) +
  geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85))# + theme(text = element_text(family="Arial"))

# linear distance traveled####

# plot with raw means:

# td<-with(behave, aggregate((total.dist.moved), list(Treatment=Treatment), mean))
# td$se<-with(behave, aggregate((total.dist.moved), list(Treatment=Treatment), 
#                               function(x) sd(x)/sqrt(length(x))))[,2]
# 
# td.plot<- ggplot(td, aes(x=Treatment, y=x, fill=Treatment)) +
#   geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
#   scale_x_discrete(limits=c("Low","Medium","High"))+
#   theme_classic() + 
#   labs(x="Risk Treatment", y="Linear Distance Traveled (mm)") +
#   theme(legend.position="none") + 
#   scale_fill_manual(values=c("grey", "grey", "grey")) +
#   theme(axis.text.x=element_text(size=20, colour="black"),
#         axis.text.y=element_text(size=20, colour="black"), 
#         axis.title=element_text(size=20))+
#   theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
#         axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
#         axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
#   theme(legend.text=element_text(size=18)) +
#   theme(legend.title =element_text(size=20))+
#   theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,205))
# td.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
#                            position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

# using plot with adjusted means instead:

dm.1e.emm.df$trt_o<-ordered(dm.1e.emm.df$Treatment,levels=c("Low","Medium","High"))

l<- ggplot(dm.1e.emm.df, aes(x=trt_o, y=emmean, fill=trt_o)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment", y ="Linear distance traveled (mm)\n(adjusted)") +
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
  scale_y_continuous(expand = c(0,0), limits = c(0, 205)) +
  geom_linerange(aes(ymin=emmean-SE, ymax=emmean+SE), size=0.5,   
                 position=position_dodge(.85))

# foraging rate (bites per minute)####

fr<-with(behave, aggregate((bites.min), list(Treatment=Treatment), mean))
fr$se<-with(behave, aggregate((bites.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

f<- ggplot(fr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment",y=(expression(atop("Foraging rate", 
                                             paste((bites~min^-1))))))+
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
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.01)) +
  geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85))# + theme(text = element_text(family="Arial"))

# interactions with conspecifics (displays per minute)####

#movements per minute with rate in parentheses

cr<-with(behave, aggregate((courtship.min), list(Treatment=Treatment), mean))
cr$se<-with(behave, aggregate((courtship.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

c <- ggplot(cr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment",y=(expression(atop("Interactions with conspecifics", 
                                             paste((displays~min^-1))))))+
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
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.151))+
  geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85))# + theme(text = element_text(family="Arial"))

# movement rate (movements per minute)####

mm<-with(behave, aggregate((movements.min), list(Treatment=Treatment), mean))
mm$se<-with(behave, aggregate((movements.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

m <- ggplot(mm, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() +
  labs(x="Risk treatment",y=(expression(Movements~min^{" -1"}))) +
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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.51),
                                                             labels = scales::number_format(accuracy = 0.01))+ #changed to 2 decimal places for movement rate)
 geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85))# + theme(text = element_text(family="Arial"))

m

# survivorship (plot only, see script for recollections for analyses) #####

# carried over from other script for recollections because I want to include all these plots in a single figure

# plotting (all trials, but not caging effects in trial 6)

# import dataset
reco<-read.csv("../2018 recollections/Data/2019.10.8.recollection.data.csv")

# data manipulation 

#adding column for survivorship, dividing recollections by 20 (initial number of fish on reefs)
reco$Survivorship<-reco$Count/20

# plot

survival<-with(reco, aggregate((Survivorship), list(Treatment=Treatment), mean))
survival$se<-with(reco, aggregate((Survivorship), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]

r <- ggplot(survival, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment", y="Survivorship") +
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
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.31),
                                                             labels = scales::number_format(accuracy = 0.01)) +
  geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85))# + theme(text = element_text(family="Arial"))

# arranging all plots together in a single figure #####

# referencing the names of the plots to be displayed on the panels

# f, r, m, c, l, e

# panels occur in this order:

# 1 2
# 3 4
# 5 6

# I've already removed x-axis labels

fr <- ggarrange(r, f, e, c, l, m, # adding null plot in the middle allows you to specify spacing between panels
                #labels = c("A)", "D)", "B)", "E)", "C)", "F)"),
                ncol = 2, nrow = 3, # widths = c(1.5, -0.09, 1.5), # changes the width of figures
                align = c("hv"),
                font.label = list(size = 16),
                hjust = -8, vjust = 0.3) # aligns labels vertically and horizontally


png("Output/figure_test_HD.png", width = 12, height = 13, units = 'in', res = 1000)
fr
dev.off()
