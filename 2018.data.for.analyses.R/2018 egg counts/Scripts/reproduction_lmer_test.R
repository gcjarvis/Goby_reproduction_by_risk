# --------------------------------------
# Description: R script for reproduction for Jarvis and Steele
# Author: George C. Jarvis
# Date: Mon Aug 10 17:17:13 2020
# Notes: what I did for egg counts (original counts): I used lmer + log likelihoods to do model selection
# then I analyzed the final model with nlme (made the most sense for df for output); so the tables in the
# paper reflect the analyses and packages that they were done with. Will likely just say that I used nlme 
# in the paper
# I'll do the same for the analyses for biomass analyses as well
# --------------------------------------

library(car)
library(tidyverse)
library(lme4)
library(nlme)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(broom.mixed)
library(emmeans) #for generating least-squares adjusted means from models (will likely help when I have to
## back-transform means if data transformation is needed)
library(sciplot)

#importing dataset####
repro<-read.csv("Data/8_10_20_big_dataset_repro.csv") #includes biomass, densities, and egg counts
head(repro)

#data for t6 t-test
t6_t_test<-read.csv("Data/t6_t_test_eggs.csv")

#data manipulation####

#adding column for average number of inhabitants throughout the trial, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$recollection_male_and_female+20)/2))

#making Year and Trial factors
repro$Year<- as.factor(repro$Year)
repro$Trial<-as.factor(repro$Trial)

#sqrt-transforming egg counts to better satisfy assumptions
repro$sqrt.egg.week<-sqrt(repro$egg.week)

#subsetting data from Trial 6 for comparison of caged vs. uncaged treatments
repro.t6<-repro[repro$Trial==6,]
#subset treatments, caged and uncaged only to test for cage effects
## Note that Treatment factor for trial 6 models should be replaced with "T6_comparison" for analyses and plotting
repro.t6<-repro.t6[repro.t6$T6_comparison == "High" | repro.t6$T6_comparison == "Uncaged", ]

#consider deleting this chunk of code
#exporting wrangled data for those that are not working with these data in R
#Data from Trials 1-5
#write.csv(repro,"Data\\egg_counts_after_data_wrangling.2020.8.10.csv", row.names = FALSE)
#Data from Trial 6 only
#write.csv(repro.t6,"Data\\Cage_effects_egg_counts_after_data_wrangling.2020.4.7.csv", row.names = FALSE)

#datasets for biomass analyses: 
#all trials, removing T-6 comparison column and then removing NA's
biomass<-select(repro, -T6_comparison)
#View(biomass)
biomass_all_cc<-drop_na(biomass)

#same but for t6 only:
biomass_t6_cc<-drop_na(repro.t6) #looks like I recollected fems from all reefs

#analyses####

#[Egg counts]

# analyzing mixed models with log-likelihood estimates and chi-square tests.
# start with full model, then reduce model first by non-significant random effects, then by NS. fixed effects
# want to definitely leave in Trial term as random effect, and Treatement, Year, T x Y, and avg.inhab as fixed effects.
# remove all NS. interactions with the covariate (fixed and random)

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#model structure for random effects: all slopes are the same (effcets are the same among trials), but intercepts are random (magnitude differs)
## Note: use "REML=F" (maximum likelihood estimates) to compare models with log-likelihood estimates (Pinheiro & Bates, 2000; Bolker et al., 2009)

# 1a. all trials, raw data for reproduction per week####
# full model
me<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me)
broom.mixed::tidy(me) #interesting code, probably won't need it here, but would come in handy 

boxplot(egg.week~Year,repro)

me$coeff

ranef(me)$Trial

coef(me)

rand(me)

levels(repro$Treatment.fact)

#same model, but with separate slopes for each trial ((avg.inhab|Treatment:Trial))
#doesn't appear to make much of a difference, and the two models are not 
#mer<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
#           (1|Trial:avg.inhab)+(avg.inhab|Treatment:Trial), REML=F, repro)
#hist(resid(mer))
#qqnorm(resid(mer))
#qqline(resid(mer))
#plot(mer)
#summary(mer)
#anova(mer)

#anova(me,mer) #not going to use this model in any of the the further analyses


#removing three-way interaction of random effect
me2<-update(me, .~. -(1|avg.inhab:Treatment:Trial))

summary(me2)
anova(me2)
Anova(me2, type = "III") # interesting that this does wald type III analysis of deviance tests, and that there are qualitative differences between this and
# - the "anova(mod...)" results

anova(me,me2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

me3<-update(me2, .~. -(1|Treatment:Trial))
summary(me3)
anova(me3)
Anova(me3, type = "III")

anova(me2,me3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

me4<-update(me3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(me4) # trial accounts for 6% of the residual variance
anova(me4)

anova(me3,me4) # P = 0.3939, doesn't suggest any difference due to the dropping of the random effect of Trial:avg.inhab

#trial effect
me5<-update(me2,.~. -(1|Trial))
summary(me5)
anova(me5)

anova(me2,me5) # no sig. effect of trial

# decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant random effects at P>0.25
## as tested for with log-likelihood tests

# me4 is the model where all random effects aside from 1|Trial have been removed

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

me6<-update(me4, .~. -(Treatment:Year:avg.inhab))
summary(me6) # trial accounts for 6% of the residual variance
anova(me6)

anova(me4,me6)# P = 0.5013

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
me7<-update(me6, .~. -(Year:avg.inhab))
summary(me7) # trial variance is still 6% of residual variance for random effects
anova(me7)

anova(me7,me6)# P = 0.779

#now want to take out Treatment:avg.inhab

me8<-update(me7, .~. -(Treatment:avg.inhab)) 
summary(me8) # trial variance is still 6% of residual variance for random effects
anova(me8)

anova(me7,me8)# P = 0.7064

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(me8))#not bad
#LS-adjusted means from model
emmeans(me8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(egg.week~Treatment,data=repro)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
me9<-update(me8, .~. -(avg.inhab))
summary(me9) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(me9)

anova(me8,me9)#chisq = 38.021      df = 1      p = 7e-10, super significant, so don't want to remove that term
#also going to include the estimate for avg.inhab from the previous model (m8)

# - Treatment*Year.fact

me10<-update(me8, .~. -(Treatment:Year))
summary(me10)
anova(me10)

anova(me8,me10) #chisq = 0.5296     df= 2     p = 0.7674

# - Treatment

me11<-update(me10,.~. -(Treatment)) #was getting 0 df when dropping treatment, 
# and it is because Treatment:Year was still in the model
summary(me11)
anova(me11)

anova(me10,me11) #chisq = 2.26     df = 2         p = 0.32 #this is acceptable, amanda has this in her paper for a factor

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

me12<-update(me10,.~. -(Year))
summary(me12)
anova(me12)

anova(me10,me12)

#using nlme, seems to make more sense with the year factor?#####
#snytax for nlme makes it hard to specify multiple random effects
#trying one different model just in case, but doesn't seem to fit my needs

#template (https://stats.stackexchange.com/questions/58669/specifying-multiple-separate-random-effects-in-lme):
#fit <- lme(Y ~ time, random=list(year=~1, date=~time), data=X, weights=varIdent(form=~1|year))

#fit <- lme(egg.week~Treatment*Year*avg.inhab, 
#           random=list(Trial=~1, Treatment=~Trial, data=repro, weights=varIdent(form=~1|Trial)))
#full model:
mod_full_lme<-lme(egg.week~Treatment*Year*avg.inhab, random = ~1|Trial, data = repro)
Anova(mod_full_lme, type = "III")
anova(mod_full_lme)
summary(mod_full_lme)
ranef(mod_full_lme)

#reduced model:
mod_final<-lme(egg.week~Treatment*Year+avg.inhab, random = ~1|Trial, data = repro)
anova(mod_final)
summary(mod_final)
ranef(mod_final)

# 1c. same analysis with log-likelihood estimates, but with sqrt.eggs as response####

#nlme:
#full model:
smod_full_lme<-lme(sqrt.egg.week~Treatment*Year*avg.inhab, random = ~1|Trial, data = repro)
anova(smod_full_lme)
summary(smod_full_lme)
ranef(smod_full_lme)
rand(me)

#reduced model:
smod_final<-lme(sqrt.egg.week~Treatment*Year+avg.inhab, random = ~1|Trial, data = repro)
Anova(smod_final, type = "III")
anova(smod_final)
summary(smod_final)
ranef(smod_final)

# removing treatment*year term (keeping this in the final model, interested in that term because it corresponds to the repro * methodology question)
smod_final1<-lme(sqrt.egg.week~Treatment+Year+avg.inhab, random = ~1|Trial, data = repro)
Anova(smod_final1, type = "III")
anova(smod_final1)
summary(smod_final1)
ranef(smod_final1)

#tried extracting MS residuals directly from the model, but I don't think this is what I want
r.rmse <- sqrt(mean(residuals(smod_final)^2))
r.rmse

#trying a different model, uncorrelated random intercept and random intercept within group
smod_full_lme<-lme(sqrt.egg.week~Treatment*Year*avg.inhab, random = ~Treatment|Trial, data = repro)
anova(smod_full_lme)
summary(smod_full_lme)
ranef(smod_full_lme)

# this is where I updated the models for reproductive output on June 13, 2021 #######

## all I have to do is change the original model
mes<-lmer(sqrt.egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
            (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, repro)
hist(resid(mes))
qqnorm(resid(mes))
qqline(resid(mes))
plot(mes)
summary(mes)
anova(mes)
rand(mes)# trial is the only random effect that seems to explain any of the variance 

# result (suggests model won't converge):

# > rand(mes)# trial is the only random effect that seems to explain any of the variance 
# boundary (singular) fit: see ?isSingular
# boundary (singular) fit: see ?isSingular
# boundary (singular) fit: see ?isSingular
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   sqrt.egg.week ~ Treatment + Year + avg.inhab + (1 | Trial) + 
#   (1 | Treatment:Trial) + (1 | Trial:avg.inhab) + (1 | avg.inhab:Treatment:Trial) + 
#   Treatment:Year + Treatment:avg.inhab + Year:avg.inhab + Treatment:Year:avg.inhab
#                                 npar  logLik    AIC     LRT Df Pr(>Chisq)
# <none>                            17 -481.44 996.89                      
# (1 | Trial)                       16 -481.76 995.53 0.63674  1     0.4249
# (1 | Treatment:Trial)             16 -481.44 994.89 0.00000  1     1.0000
# (1 | Trial:avg.inhab)             16 -481.44 994.89 0.00000  1     0.9990
# (1 | avg.inhab:Treatment:Trial)   16 -481.44 994.89 0.00000  1     1.0000

Anova(mes, type = "III")

# result:

# > Anova(mes, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
# Chisq Df Pr(>Chisq)    
# (Intercept)               5.5189  1    0.01881 *  
# Treatment                 0.3339  2    0.84626    
# Year                      0.0211  1    0.88446    
# avg.inhab                27.6919  1  1.423e-07 ***
# Treatment:Year            1.9908  2    0.36958    
# Treatment:avg.inhab       0.6463  2    0.72385    
# Year:avg.inhab            0.0000  1    0.99878    
# Treatment:Year:avg.inhab  1.9460  2    0.37795    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# also, model won't even converge if I keep all of these higher-level interactions in it

# proper model reduction, after removing the nonsignificant interactions with the random effects

mes1<-lmer(sqrt.egg.week ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) +
             (Year*avg.inhab) + (Treatment*Year*avg.inhab) + (1|Trial), REML=F, repro)
Anova(mes1, type = "III")


# result:

# > Anova(mes1, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
#                               Chisq Df Pr(>Chisq)    
#   (Intercept)               5.5189  1    0.01881 *  
#   Treatment                 0.3339  2    0.84626    
#   Year                      0.0211  1    0.88446    
#   avg.inhab                27.6918  1  1.423e-07 ***
#   Treatment:Year            1.9908  2    0.36958    
#   Treatment:avg.inhab       0.6463  2    0.72385    
#   Year:avg.inhab            0.0000  1    0.99878    
#   Treatment:Year:avg.inhab  1.9460  2    0.37795    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# removing trial, running log-likelihood test

mes1a<-lm(sqrt.egg.week ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) +
             (Year*avg.inhab) + (Treatment*Year*avg.inhab), repro)
Anova(mes1a, type = "III")
summary(mes1a)

2*(logLik(mes1) - logLik(mes1a)) # Chi2 =  0.67
pchisq(2*(logLik(mes1) - logLik(mes1a)), df = 1, lower.tail=F) # P = 0.41, for df = 1

# reporting that for effect of trial in model, moving on to rest of model reduction

# removing 3-way interaction:

mes2<-lmer(sqrt.egg.week ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) +
             (Year*avg.inhab) + (1|Trial), REML=F, repro)
Anova(mes2, type = "III")

# result:

# > Anova(mes2, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
# Chisq Df Pr(>Chisq)    
# (Intercept)          4.4915  1    0.03406 *  
#   Treatment            0.1092  2    0.94686    
# Year                 0.0913  1    0.76257    
# avg.inhab           25.6124  1  4.173e-07 ***
#   Treatment:Year       0.1182  2    0.94261    
# Treatment:avg.inhab  0.2009  2    0.90445    
# Year:avg.inhab       0.2144  1    0.64334    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# testing chi-square values with log-likelihood tests, thought df = 1, but df's here are 1 or 2 (for multi-level factors)

2*(logLik(mes1) - logLik(mes2))# Chi2 =  1.928381, very similar to what I got for the three-way interaction term
pchisq(2*(logLik(mes1) - logLik(mes2)), df = 1, lower.tail=F) # P = 0.165, which is different from what you get if you use df = 2
pchisq(2*(logLik(mes1) - logLik(mes2)), df = 2, lower.tail=F) # P = 0.381, which is more similar to what the output table looks like

# just going to report the analysis of deviance from the car package...moving on with model selection

# removing year*inhab

mes3<-lmer(sqrt.egg.week ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) +
              (1|Trial), REML=F, repro)
Anova(mes3, type = "III")

# result:

# > Anova(mes3, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
#                         Chisq Df Pr(>Chisq)    
#   (Intercept)          6.3916  1    0.01147 *  
#   Treatment            0.1143  2    0.94446    
#   Year                 1.0521  1    0.30503    
#   avg.inhab           33.7288  1  6.335e-09 ***
#   Treatment:Year       0.1275  2    0.93826    
#   Treatment:avg.inhab  0.2238  2    0.89415    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# removing treatment*inhab

mes4<-lmer(sqrt.egg.week ~ Treatment + Year + avg.inhab + (Treatment*Year) +
             (1|Trial), REML=F, repro)
Anova(mes4, type = "III")

# results:

# > Anova(mes4, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
# Chisq Df Pr(>Chisq)    
# (Intercept)     6.9276  1   0.008487 ** 
#   Treatment       2.2648  2   0.322261    
# Year            1.1257  1   0.288695    
# avg.inhab      36.7713  1  1.328e-09 ***
#   Treatment:Year  0.1352  2   0.934647    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# removing treatment*year (worth reporting this NS result in the paper to say that tagging method didn't ultimately end up affecting reproduction...
# but that I would highly recommend not anesthitizing in MS222)

# removing treatment*year

mes5<-lmer(sqrt.egg.week ~ Treatment + Year + avg.inhab +
             (1|Trial), REML=F, repro)
Anova(mes5, type = "III")

# results:

# > Anova(mes5, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
# Chisq Df Pr(>Chisq)    
# (Intercept)  6.8622  1   0.008804 ** 
#   Treatment    2.2600  2   0.323034    
# Year         1.1866  1   0.276018    
# avg.inhab   36.9198  1  1.231e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# let's see what happens if I remove the year effect (again, worth reporting that there was no statistical difference between the two years)

mes6<-lmer(sqrt.egg.week ~ Treatment + avg.inhab +
             (1|Trial), REML=F, repro)
Anova(mes6, type = "III")

# putting random interactions back in to see if this runs now (still has convergence issues, even whe I reduce all the other random interactions)
# 
# mes6<-lmer(sqrt.egg.week ~ Treatment + avg.inhab +
#              (1|Trial) + (1|Treatment:Trial) +
#              (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, repro)
# Anova(mes6, type = "III")
# summary(mes6)


# result:

# > Anova(mes6, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
#               Chisq Df Pr(>Chisq)    
# (Intercept)  6.4865  1    0.01087 *  
# Treatment    2.3141  2    0.31441    
# avg.inhab   35.7788  1   2.21e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# coefficient for avg.inhab?

summary(mes6)

# result:

> summary(mes6)
# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: sqrt.egg.week ~ Treatment + avg.inhab + (1 | Trial)
#    Data: repro
# 
#      AIC      BIC   logLik deviance df.resid 
#    978.5    994.7   -483.2    966.5      104 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.3520 -0.5660  0.1177  0.6636  2.3516 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Trial    (Intercept)  18.48    4.299  
#  Residual             369.77   19.229  
# Number of obs: 110, groups:  Trial, 6
# 
# Fixed effects:
#             Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)  -36.146     14.192  96.538  -2.547   0.0125 *  
# Treatment1     1.101      2.554 105.192   0.431   0.6673    
# Treatment2     2.807      2.629 103.981   1.067   0.2882    
# avg.inhab      6.552      1.095 104.186   5.982 3.15e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#            (Intr) Trtmn1 Trtmn2
# Treatment1 -0.041              
# Treatment2 -0.015 -0.484       
# avg.inhab  -0.984  0.034  0.019


AIC(mes5, mes6)

# result:

# > AIC(mes5, mes6)
# df      AIC
# mes5  7 979.3857
# mes6  6 978.4844

# model without year is slightly better than the model with year in it, but not by much

# okay, let's look at the emmeans (and plot?) for the fully-reduced model (mes6)

emmeans(mes6, pairwise~Treatment)# seems like this is the best one to go with (got an error message when I tried to include year and avg.inhab)

# result:

# > emmeans(mes6, pairwise~Treatment)
# $emmeans
# Treatment emmean   SE   df lower.CL upper.CL
# High        48.4 3.77 25.6     40.6     56.2
# Low         50.1 3.93 30.7     42.1     58.1
# Medium      43.4 3.93 30.8     35.4     51.4
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate   SE  df t.ratio p.value
# High - Low       -1.71 4.54 108 -0.376  0.9252 
# High - Medium     5.01 4.54 108  1.102  0.5149 
# Low - Medium      6.71 4.66 107  1.439  0.3245 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# based on LS- means: low > high > med

# for estimates and std. error for plots: (mean - se, mean, mean + se)^2 (each term +/- sem is squared, because these are square root-transformed means)

# high-risk

# plot for this:

mes6.emm <- emmeans(mes6, ~ Treatment)
plot(mes6.emm)
plot(warp.emm, by = NULL, comparisons = TRUE, adjust = "mvt", 
     horizontal = FALSE, colors = "darkgreen")


# how does this compare to plots of the raw data?

bargraph.CI(x.factor = Treatment, response = sqrt.egg.week, main="raw data - sqrt repro/reef/week x risk treatment", 
            xlab="Treatment", ylab="reproduction per reef per week", data = repro)

# looks the same as the plots for emmeans, how about the plots with raw data for egg counts?

bargraph.CI(x.factor = Treatment, response = egg.week, main="raw data - egg.week (not sqrt)/reef/week x risk treatment", 
            xlab="Treatment", ylab="reproduction per reef per week", data = repro)

# decided to go with data (means +/- sem) from emmeans

# exporting results as a df to make back transforming easier

egg.emm.df <- as.data.frame(egg.emm$emmeans)
egg.emm.df

# results (seems like more precise values for means and se):

# > egg.emm.df
# Treatment   emmean       SE       df lower.CL upper.CL
# 1      High 48.39883 3.769955 25.61806 40.64395 56.15371
# 2       Low 50.10488 3.927950 30.69713 42.09057 58.11919
# 3    Medium 43.39056 3.930325 30.81328 35.37264 51.40848

# adding columns for back-transformed means, + se, and - se

egg.emm.df <- egg.emm.df %>% mutate(mean = emmean^2, up_se = (((emmean + SE)^2) - mean), low_se = (mean - ((emmean - SE)^2)))

egg.emm.df

# table with means and se back-transformed:

# > egg.emm.df
# Treatment   emmean       SE       df lower.CL upper.CL     mean    up_se   low_se
# 1      High 48.39883 3.769955 25.61806 40.64395 56.15371 2342.447 379.1353 350.7102
# 2       Low 50.10488 3.927950 30.69713 42.09057 58.11919 2510.499 409.0477 378.1901
# 3    Medium 43.39056 3.930325 30.81328 35.37264 51.40848 1882.741 356.5254 325.6305

# plotting repro per reef per week

egg.emm.df$trt_o<-ordered(egg.emm.df$Treatment,levels=c("Low","Medium","High"))

repro_adjusted_plot <- ggplot(egg.emm.df, aes(x=trt_o, y=mean, fill=trt_o)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk treatment", y =expression(atop(Reproduction~reef^{" -1"}~week^{"-1"},paste("(adjusted)")))) +
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
  scale_y_continuous(expand = c(0,0), limits = c(0, 3050)) +
  geom_linerange(aes(ymin=mean-low_se, ymax=mean+up_se), size=0.5,   
                           position=position_dodge(.85))


# exporting figure for reproduction

pdf("Output/repro_ms_fig_barplot_adjusted_means_for_density.pdf", width = 8, height = 7)

repro_adjusted_plot

dev.off()



# let's look at the emmeans for the model with year and avg. inhab in it

emmeans(mes5, pairwise~Treatment)

# results (seems like convergence issue?):

# > emmeans(mes5, pairwise~Treatment)
# $emmeans
# Treatment emmean   SE   df lower.CL upper.CL
# High        48.2 3.84 37.7     40.4     55.9
# Low         50.2 3.99 44.8     42.1     58.2
# Medium      43.5 4.00 45.0     35.4     51.5
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate   SE  df t.ratio p.value
# High - Low       -2.00 4.54 108 -0.441  0.8984 
# High - Medium     4.72 4.54 108  1.039  0.5543 
# Low - Medium      6.72 4.66 107  1.442  0.3233 
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 
# 
# Warning message:
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00498512 (tol = 0.002, component 1)

#reduced model:
mes_red<-lmer(sqrt.egg.week ~ Treatment*Year+avg.inhab + (1|Trial), REML=F, repro)
hist(resid(mes_red))
qqnorm(resid(mes_red))
qqline(resid(mes_red))
plot(mes_red)
summary(mes_red)
rand(mes_red) # can report the log likelihood value for the random effect of trial

# result:

# > rand(mes_red) # can report the log likelihood value for the random effect of trial
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   sqrt.egg.week ~ Treatment + Year + avg.inhab + (1 | Trial) + 
#   Treatment:Year
# npar  logLik    AIC     LRT Df Pr(>Chisq)
# <none>         9 -482.63 983.25                      
# (1 | Trial)    8 -482.99 981.97 0.72316  1     0.3951

anova(mes_red)

# result:

# > anova(mes_red)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# Treatment        835.1   417.6     2 104.221  1.1324    0.3262    
# Year             415.1   415.1     1   5.608  1.1257    0.3322    
# avg.inhab      13559.1 13559.1     1 101.954 36.7713 2.269e-08 ***
# Treatment:Year    49.8    24.9     2 103.984  0.0676    0.9347    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now with car package for good measure

Anova(mes_red, type = "III")

# result:

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
# Chisq Df Pr(>Chisq)    
# (Intercept)     6.9276  1   0.008487 ** 
#   Treatment       2.2648  2   0.322261    
# Year            1.1257  1   0.288695    
# avg.inhab      36.7713  1  1.328e-09 ***
#   Treatment:Year  0.1352  2   0.934647    

# same qualitative result, no treatment x year effect, so removing = the effect of treatment on reproduction did not depend on the tagging procedure we used

mes_red__no_int<-lmer(sqrt.egg.week ~ Treatment + Year + avg.inhab + (1|Trial), REML=F, repro)
summary(mes_red_no_int)





#note: this is only the case with sqrt-transformed data

anova(me,mes)# seems like the log-transformed data might be a better model, based on lower AIC

#removing three-way interaction of random effect
mes2<-update(mes, .~. -(1|avg.inhab:Treatment:Trial))
summary(mes2)
anova(mes2)

anova(mes,mes2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

mes3<-update(mes2, .~. -(1|Treatment:Trial))
summary(mes3)
anova(mes3)

anova(mes2,mes3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

mes4<-update(mes3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(mes4) # trial accounts for 4% of the residual variance
anova(mes4)

anova(mes3,mes4) # P = 1, which is different from the analyses for raw data

#trial effect
mes5<-update(mes2,.~. -(1|Trial))
summary(mes5)
anova(mes5)

anova(mes2,mes5) # no sig. effect of trial

# decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant random effects at P>0.25
## as tested for with log-likelihood tests

# me4 is the model where all random effects aside from 1|Trial have been removed

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

mes6<-update(mes4, .~. -(Treatment:Year:avg.inhab))
summary(mes6) # trial accounts for 4% of the residual variance
anova(mes6)

anova(mes4,mes6)# P = 0.42

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
mes7<-update(mes6, .~. -(Year:avg.inhab))
summary(mes7) # trial variance is still 6% of residual variance for random effects
anova(mes7)

anova(mes7,mes6)# P = 0.91

#now want to take out Treatment:avg.inhab

mes8<-update(mes7, .~. -(Treatment:avg.inhab)) 
summary(mes8) # trial variance is still 6% of residual variance for random effects
anova(mes8)

anova(mes7,mes8)# P = 0.86

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(mes8))#not bad
#LS-adjusted means from model
emmeans(mes8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(sqrt.egg.week~Treatment,data=repro)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
mes9<-update(mes8, .~. -(avg.inhab))
summary(mes9) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(mes9)

anova(mes8,mes9)#chisq = 27.78      df = 1      p = 1.36e-10, super significant, so don't want to remove that term
#also going to include the estimate for avg.inhab from the previous model (mes8)

# after doing log-likelihood tests, only random effect that explained any variance was 
# trial, so I'm going to rerun the reduced model with that, and the non-significant
# interactions between fixed effects



# - Treatment*Year.fact

mes10<-update(mes8, .~. -(Treatment:Year))
summary(mes10)
anova(mes10)

anova(mes8,mes10) #chisq = 0.193     df= 2     p = 0.908

# - Treatment

mes11<-update(mes10,.~. -(Treatment)) #was getting 0 df when dropping treatment, 
# and it is because Treatment:Year was still in the model
summary(mes11)
anova(mes11)

anova(mes10,mes11) #chisq = 1.79    df = 2         p = 0.41 #this is acceptable, amanda has this in her paper for a factor

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

mes12<-update(mes10,.~. -(Year))
summary(mes12)
anova(mes12)

anova(mes10,mes12) #chisq = 0.002    df = 2         p = 0.41

# 8-22-2020: sqrt-data with nlme model####

#trial as fixed facor
tf<-aov(sqrt.egg.week~Treatment*Year*avg.inhab*Trial, data=repro)
hist(resid(tf))
plot(tf)
anova(tf)
summary(tf)

me1<-update(tf, .~. -(Treatment:avg.inhab:Trial))
summary(me1)
anova(me1)

me2<-update(me1, .~. -(Treatment:Year:avg.inhab:Trial))
anova(me2)

me3<-update(me2, .~. -(Treatment:Year:avg.inhab))
anova(me3)

me4<-update(me3, .~. -(avg.inhab:Trial))
anova(me4)

me5<-update(me4, .~. -(Year:avg.inhab:Trial))
anova(me5)

me6<-update(me5, .~. -(Treatment:Trial))
anova(me6)

me7<-update(me6, .~. -(Treatment:Year:Trial))
anova(me7)

me8<-update(me7, .~. -(Year:avg.inhab))
anova(me8)

me9<-update(me8, .~. -(Treatment:avg.inhab))
anova(me9)

boxplot(sqrt.egg.week~Trial, data=repro)

me6<-update(me5, .~. -(Treatment:Year))
anova(me6)


#using nlme, seems to make more sense with the year factor?#####
sqmod_final<-lme(sqrt.egg.week~Treatment*Year+avg.inhab, random = ~1 + Year|Trial, data = repro)
anova(sqmod_final)
summary(sqmod_final)
ranef(sqmod_final)

# 1d. Trial 6 data only, comparing square root data for high-risk caged to uncaged treatments in Trial 6 with ANCOVA####
#note: Year is not needed anymore, and Trial is removed as well

View(repro.t6)

mes13<-lm(sqrt.egg.week~T6_comparison*avg.inhab,data=repro.t6)
hist(resid(mes13)) # not so great, but might look better after data transformation
qqnorm(resid(mes13))
qqline(resid(mes13))

anova(mes13)

# removing non-significant interaction of covariate

mes14<-lm(sqrt.egg.week~T6_comparison+avg.inhab,data=repro.t6)
hist(resid(mes14)) #not a super great fit, but I think it's okay
qqnorm(resid(mes14))
qqline(resid(mes14))

anova(mes13,mes14) # no difference in model after interaction was removed

anova(mes14)

#removing covariate
mes15<-lm(sqrt.egg.week~T6_comparison,data=repro.t6)
hist(resid(mes15)) #not a super great fit, but I think it's okay
qqnorm(resid(mes15))
qqline(resid(mes15))

anova(mes15)

#t test for trial 6
t.test(t6_t_test$sqrt_high_cage ,t6_t_test$sqrt_high_uncaged)

#getting means and se for trial 6 only

#generic code using df
t6_means<-with(repro.t6, aggregate((sqrt.egg.week), list(T6_comparison=T6_comparison), mean))
t6_means
#now apply the se function to the 4th column [,3]
t6_means$se<-with(repro.t6, aggregate((sqrt.egg.week), list(T6_comparison=T6_comparison), function(x) sd(x)/sqrt(length(x))))[,2]
t6_means

#not sure how to back-transform standard error, but can back-transform the means to report in the paper
t6_means$btm<-(t6_means$x)^2
t6_means

boxplot(egg.week~T6_comparison,data = repro.t6)

(-0.56836)^2

?qf




#bimodal dist, tried a tranformation, but I don't think this is legit, and I don't have a good justification for doing it
transformed <- abs(repro.t6$sqrt.egg.week - mean(repro.t6$sqrt.egg.week))
mes15<-lm(transformed~T6_comparison,data=repro.t6)
hist(resid(mes15)) #dist does look better though
qqnorm(resid(mes15))
qqline(resid(mes15))
anova(mes15)

anova(mes15,mes14)


# LS-means for treatment from model:
emmeans(mes14, pairwise~T6.comparison)
boxplot(sqrt.egg.week~T6.comparison,data=repro.t6)# variances don't look too bad

#plotting ANCOVA figure for RAW data; this is what I present in the manuscript####

#Order the treatments: Low --> Med --> High
repro$Treatment.ord<-ordered(repro$Treatment,levels=c("Low","Medium","High"))

#png("Output/2020.4.11.gobies.per.reef.9.5x5.5.300dpi.png", width = 9.5, height = 5.7, units = 'in', res = 300)

anc1<-ggplot(repro, aes(avg.inhab, egg.week, shape=Treatment.ord, linetype=Treatment.ord, col=Treatment.ord)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x=(bquote('Gobies '~Reef^ -1*'')),y=(bquote('Reproducton '~Reef^ -1~ Week^-1*'')))+
  expand_limits(y=0)+
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=16)) +
  theme(legend.title =element_text(size=17))+
  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,12000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
anc1

#dev.off()

# rerunning same model, but with trial as fixed factor ####

fixed_T<-aov(egg.week ~ Treatment*Year*avg.inhab*(Trial%in%Year), repro)
hist(resid(fixed_T))
qqnorm(resid(fixed_T))
qqline(resid(fixed_T))
plot(fixed_T)
summary(fixed_T)
anova(fixed_T)

boxplot(egg.week~Treatment,repro)

me2<-update(fixed_T, .~. -(Treatment:Year:avg.inhab:Trial))
summary(me2)
anova(me2)

anova(fixed_T,me2)

me3<-update(me3, .~. -(Year:avg.inhab:Trial))
summary(me3)
anova(me3)

anova(me2,me3) #error?

me4<-update(me3, .~. -(Treatment:Year:avg.inhab))
summary(me4)
anova(me4)

anova(me3,me4)

me5<-update(me4, .~. -(Year:avg.inhab))
summary(me5)
anova(me5)

anova(me4,me5)

m1<-aov(egg.week ~ Treatment*Year+avg.inhab+(Trial%in%Year), repro)
hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
plot(m1)
summary(m1)
anova(m1)


m1a<-aov(egg.week ~ Treatment*Year+avg.inhab+(Trial%in%Year)+Treatment*(Trial%in%Year), repro)
hist(resid(m1a))
qqnorm(resid(m1a))
qqline(resid(m1a))
plot(m1a)
summary(m1a)
anova(m1a)

anova(m1,m1a)

boxplot(egg.week~(Trial%in%Year), repro)

boxplot(egg.week~Treatment*(Trial%in%Year), repro)

# runnning per capita output as response variable with same model
head(repro1)

mcap<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female+(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

boxplot(reproduction_per_capita_female_biomass_per_day_per_gram~(Trial%in%Year), repro1)

mcap1<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female+Treatment*(Trial%in%Year), repro1)
hist(resid(mcap1))
qqnorm(resid(mcap1))
qqline(resid(mcap1))
plot(mcap1)
summary(mcap1)
anova(mcap1) #no trial x treatment effect

boxplot(recollection_female~Treatment, repro1)

mcap2<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female
           +Treatment*(Trial%in%Year)+Treatment*recollection_female, repro1)
hist(resid(mcap2))
qqnorm(resid(mcap2))
qqline(resid(mcap2))
plot(mcap2)
summary(mcap2)
anova(mcap2)

# findings: 1) sig. differences by trial when included as a nested fixed factor. 
# 2) Same qualitative results for per capita reproduction per female biomass per day (reproduction per female per gram per day)
# 3) 

#trial as random
mcap3<-lmer(sqrt(reproduction_per_capita_female_biomass_per_day_per_gram) ~ Treatment*Year*recollection_female*(1|Trial)
            +(1|Treatment:Trial)+(1|recollection_female:Trial)+(1|recollection_female:Trial:Treatment), REML = F, repro)
hist(resid(mcap3))
qqnorm(resid(mcap3))
qqline(resid(mcap3))
plot(mcap3)
summary(mcap3)
anova(mcap3)

rand(mcap3)

#dropping all non-significant random effects
mcap4<-lmer(sqrt(reproduction_per_capita_female_biomass_per_day_per_gram) ~ Treatment*Year*recollection_female*(1|Trial), REML = F, repro)
hist(resid(mcap4))
qqnorm(resid(mcap4))
qqline(resid(mcap4))
plot(mcap4)
summary(mcap4)
anova(mcap4)

#adding covariate, not crossing
mcap5<-lmer(sqrt(reproduction_per_capita_female_biomass_per_day_per_gram) ~ Treatment*Year+recollection_female*(1|Trial), REML = F, repro)
hist(resid(mcap5))
qqnorm(resid(mcap5))
qqline(resid(mcap5))
plot(mcap5)
summary(mcap5)
anova(mcap5)

boxplot(sqrt(reproduction_per_capita_female_biomass_per_day_per_gram)~Treatment, data = repro)

### biomass recollected

mcap<-aov(per_capita_female_biomass ~ Treatment*Year+recollection_female+(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

mcap<-aov(per_capita_female_biomass ~ Treatment*Year+recollection_female+Treatment*(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

boxplot(per_capita_female_biomass~Treatment, repro1)

mcap<-lmer(per_capita_female_biomass ~ Treatment*Year+recollection_female+(1|Trial), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

# 8/2/20 rerunning models with avg_den_24h and avg_den_seen_max ####

#quick model to see if there were differences in avg_den_24hr

boxplot(avg_den_24h~Treatment, repro) # doesn't seem like it, which is interesting when
# we consider using it as a covariate. Might be an argument to remove it as a covariate?

re<-lm(avg_den_24h~Treatment*Year*Density*(Trial %in% Year), repro)
hist(resid(re))
qqnorm(resid(re))
anova(re)

me1<-update(re, .~. -(Treatment:Year:Density:Trial))
summary(me1)
anova(me1)

me2<-update(me1, .~. -(Year:Density:Trial))
anova(me2)

me3<-update(me2, .~. -(Treatment:Year:Trial))
anova(me3)

me4<-update(me3, .~. -(Treatment:Year:Density))
anova(me4)

me5<-update(me4, .~. -(Treatment:avg_den_max_seen))
anova(me5)

me6<-update(me5, .~. -(Treatment:Year))
anova(me6)

emmeans(me6, pairwise~Treatment)

me<-lmer(sqrt.egg.week ~ Treatment*Year*avg_den_24h + (1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg_den_24h)+(1|avg_den_24h:Treatment:Trial), REML=F, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me) #hmmm, looks like there's a three-way interaction now?

boxplot(egg.week~avg_den_24h,repro)

me$coeff

ranef(me)$Trial

coef(me)

rand(me)

levels(repro$Treatment)

me<-lmer(egg.week ~ Treatment*Year*avg_den_max_seen + (1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg_den_max_seen)+(1|avg_den_max_seen:Treatment:Trial), REML=F, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me) #hmmm, looks like there's a three-way interaction now?

#trying trial as fixed

me<-lm(egg.week ~ Treatment*Year*avg_den_max_seen*(Trial %in% Year), repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me) 

me<-lm(egg.week ~ Treatment*Year*avg_den_max_seen+(Trial %in% Year), repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me) 

mes<-lm(sqrt.egg.week ~ Treatment*Year*avg_den_max_seen+(Trial %in% Year), repro)
hist(resid(mes))
qqnorm(resid(mes))
qqline(resid(mes))
plot(mes)
summary(mes)
anova(mes) 

anova(me,mes)


#removing all random effects and rerunning

me<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg_den_24h)+(1|avg_den_24h:Treatment:Trial), REML=F, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me) #hmmm, looks like there's a three-way interaction now?

boxplot(egg.week~Treatment*Year,repro)

me1<-update(me, .~. -(Treatment:Year:avg_den_max_seen:Trial))
summary(me1)
anova(me1)

me2<-update(me1, .~. -(Treatment:Year:avg_den_max_seen))
anova(me2)

me3<-update(me2, .~. -(Year:Trial))
anova(me3)

me4<-update(me3, .~. -(Year:avg_den_max_seen ))
anova(me4)

me5<-update(me4, .~. -(Treatment:avg_den_max_seen))
anova(me5)

me6<-update(me5, .~. -(Treatment:Year))
anova(me6)

emmeans(me6, pairwise~Treatment)

me<-lmer(egg.week ~ Treatment*Year*avg_den_max_seen + (1|Trial), REML=F, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me)

#same model, but with separate slopes for each trial ((avg.inhab|Treatment:Trial))
#doesn't appear to make much of a difference, and the two models are not 
#mer<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
#           (1|Trial:avg.inhab)+(avg.inhab|Treatment:Trial), REML=F, repro)
#hist(resid(mer))
#qqnorm(resid(mer))
#qqline(resid(mer))
#plot(mer)
#summary(mer)
#anova(mer)

#anova(me,mer) #not going to use this model in any of the the further analyses


#removing three-way interaction of random effect
me2<-update(me, .~. -(1|avg.inhab:Treatment:Trial))
summary(me2)
anova(me2)

anova(me,me2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

me3<-update(me2, .~. -(1|Treatment:Trial))
summary(me3)
anova(me3)

anova(me2,me3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

me4<-update(me3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(me4) # trial accounts for 6% of the residual variance
anova(me4)

anova(me3,me4) # P = 0.3939, doesn't suggest any difference due to the dropping of the random effect of Trial:avg.inhab

#trial effect
me5<-update(me2,.~. -(1|Trial))
summary(me5)
anova(me5)

anova(me2,me5) # no sig. effect of trial

# decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant random effects at P>0.25
## as tested for with log-likelihood tests

# me4 is the model where all random effects aside from 1|Trial have been removed

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

me6<-update(me4, .~. -(Treatment:Year:avg.inhab))
summary(me6) # trial accounts for 6% of the residual variance
anova(me6)

anova(me4,me6)# P = 0.5013

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
me7<-update(me6, .~. -(Year:avg.inhab))
summary(me7) # trial variance is still 6% of residual variance for random effects
anova(me7)

anova(me7,me6)# P = 0.779

#now want to take out Treatment:avg.inhab

me8<-update(me7, .~. -(Treatment:avg.inhab)) 
summary(me8) # trial variance is still 6% of residual variance for random effects
anova(me8)

anova(me7,me8)# P = 0.7064

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(me8))#not bad
#LS-adjusted means from model
emmeans(me8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(egg.week~Treatment,data=repro)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
me9<-update(me8, .~. -(avg.inhab))
summary(me9) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(me9)

anova(me8,me9)#chisq = 38.021      df = 1      p = 7e-10, super significant, so don't want to remove that term
#also going to include the estimate for avg.inhab from the previous model (m8)

# - Treatment*Year.fact

me10<-update(me8, .~. -(Treatment:Year))
summary(me10)
anova(me10)

anova(me8,me10) #chisq = 0.5296     df= 2     p = 0.7674

# - Treatment

me11<-update(me10,.~. -(Treatment)) #was getting 0 df when dropping treatment, 
# and it is because Treatment:Year was still in the model
summary(me11)
anova(me11)

anova(me10,me11) #chisq = 2.26     df = 2         p = 0.32 #this is acceptable, amanda has this in her paper for a factor

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

me12<-update(me10,.~. -(Year))
summary(me12)
anova(me12)

anova(me10,me12)

# 1b. Trial 6 data only, comparing high-risk caged to uncaged treatments in Trial 6 with ANCOVA####
#note: Year is not needed anymore, and Trial is removed as well

me13<-lm(egg.week~T6_comparison*avg.inhab,data=repro.t6)
hist(resid(me13)) # not so great, but might look better after data transformation
qqnorm(resid(me13))
qqline(resid(me13))

anova(me13)

# removing non-significant interaction of covariate

me14<-lm(egg.week~T6.comparison+avg.inhab,data=repro.t6)
hist(resid(me14)) #not a super great fit, but I think it's okay
qqnorm(resid(me14))
qqline(resid(me14))

anova(me13,me14) # no difference in model after interaction was removed

# LS-means for treatment from model:
emmeans(me14, pairwise~T6.comparison)
boxplot(egg.week~T6.comparison,data=repro.t6)# variances don't look too bad

# 1c. same analysis with log-likelihood estimates, but with sqrt.eggs as response####
## all I have to do is change the orignial model
mes<-lmer(sqrt.egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
            (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, repro)
hist(resid(mes))
qqnorm(resid(mes))
qqline(resid(mes))
plot(mes)
summary(mes)
anova(mes)

rand(mes)# trial is the only random effect that seems to exaplain any of the variance

#note: this is only the case with sqrt-transformed data

anova(me,mes)# seems like the log-transformed data might be a better model, based on lower AIC

#removing three-way interaction of random effect
mes2<-update(mes, .~. -(1|avg.inhab:Treatment:Trial))
summary(mes2)
anova(mes2)

anova(mes,mes2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

mes3<-update(mes2, .~. -(1|Treatment:Trial))
summary(mes3)
anova(mes3)

anova(mes2,mes3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

mes4<-update(mes3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(mes4) # trial accounts for 4% of the residual variance
anova(mes4)

anova(mes3,mes4) # P = 1, which is different from the analyses for raw data

#trial effect
mes5<-update(mes2,.~. -(1|Trial))
summary(mes5)
anova(mes5)

anova(mes2,mes5) # no sig. effect of trial

# decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant random effects at P>0.25
## as tested for with log-likelihood tests

# me4 is the model where all random effects aside from 1|Trial have been removed

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

mes6<-update(mes4, .~. -(Treatment:Year:avg.inhab))
summary(mes6) # trial accounts for 4% of the residual variance
anova(mes6)

anova(mes4,mes6)# P = 0.42

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
mes7<-update(mes6, .~. -(Year:avg.inhab))
summary(mes7) # trial variance is still 6% of residual variance for random effects
anova(mes7)

anova(mes7,mes6)# P = 0.91

#now want to take out Treatment:avg.inhab

mes8<-update(mes7, .~. -(Treatment:avg.inhab)) 
summary(mes8) # trial variance is still 6% of residual variance for random effects
anova(mes8)

anova(mes7,mes8)# P = 0.86

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(mes8))#not bad
#LS-adjusted means from model
emmeans(mes8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(sqrt.egg.week~Treatment,data=repro)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
mes9<-update(mes8, .~. -(avg.inhab))
summary(mes9) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(mes9)

anova(mes8,mes9)#chisq = 27.78      df = 1      p = 1.36e-10, super significant, so don't want to remove that term
#also going to include the estimate for avg.inhab from the previous model (mes8)

# after doing log-likelihood tests, only random effect that explained any variance was 
# trial, so I'm going to rerun the reduced model with that, and the non-significant
# interactions between fixed effects



# - Treatment*Year.fact

mes10<-update(mes8, .~. -(Treatment:Year))
summary(mes10)
anova(mes10)

anova(mes8,mes10) #chisq = 0.193     df= 2     p = 0.908

# - Treatment

mes11<-update(mes10,.~. -(Treatment)) #was getting 0 df when dropping treatment, 
# and it is because Treatment:Year was still in the model
summary(mes11)
anova(mes11)

anova(mes10,mes11) #chisq = 1.79    df = 2         p = 0.41 #this is acceptable, amanda has this in her paper for a factor

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

mes12<-update(mes10,.~. -(Year))
summary(mes12)
anova(mes12)

anova(mes10,mes12) #chisq = 0.002    df = 2         p = 0.41

# 1d. Trial 6 data only, comparing square root data for high-risk caged to uncaged treatments in Trial 6 with ANCOVA####
#note: Year is not needed anymore, and Trial is removed as well

mes13<-lm(sqrt.egg.week~T6.comparison*avg.inhab,data=repro.t6)
hist(resid(mes13)) # not so great, but might look better after data transformation
qqnorm(resid(mes13))
qqline(resid(mes13))

anova(mes13)

# removing non-significant interaction of covariate

mes14<-lm(sqrt.egg.week~T6.comparison+avg.inhab,data=repro.t6)
hist(resid(mes14)) #not a super great fit, but I think it's okay
qqnorm(resid(mes14))
qqline(resid(mes14))

anova(mes13,mes14) # no difference in model after interaction was removed

anova(mes14)

# LS-means for treatment from model:
emmeans(mes14, pairwise~T6.comparison)
boxplot(sqrt.egg.week~T6.comparison,data=repro.t6)# variances don't look too bad

# [reproduction in response to changes in biomass] ####
# no statistical difference in the biomass of fish recollected among treatments

#1. tests for differences in biomass among treatments:
#already did the lon-likelihood tests on my Monash computer, so I know what the final model should be
# - will just use the same procedure as I did for egg counts, and use the ran(biomass_full)
head(repro)

#full modelfor log-likelihood tests:
biomass_full<-lmer(per_capita_female_biomass~Treatment*Year*recollection_female + (1|Trial)+
                   (1|Treatment:Trial) +(1|Trial:recollection_female)+
                   (1|recollection_female:Treatment:Trial),
                   REML=F, biomass_all_cc)
hist(resid(biomass_full))
rand(biomass_full)
anova(biomass_full)

boxplot(per_capita_female_biomass~Treatment*Trial,data = repro)
boxplot(recollection_female~Treatment*Trial,data = repro)

#reduced model lmer
biomass_red_lmer<-lmer(per_capita_female_biomass~Treatment*Year*recollection_female 
                       +(1|Trial:recollection_female)+ (1|Trial),
                   REML=F, biomass_all_cc)
anova(biomass_red_lmer)
summary(biomass_red_lmer)

#full model with nlme for supplementary table
#had to make a single variable for trial*number female recollected ("tnf", see below model)
biomass_full_nlme<-lme(per_capita_female_biomass~Treatment*Year*recollection_female, 
                random= ~1,biomass_all_cc, varIdent(form=~1|Trial))
hist(resid(biomass_full_nlme))
anova(biomass_full_nlme)
summary(biomass_full_nlme)

biomass_all_cc$tn<-as.numeric(biomass_all_cc$Trial)
biomass_all_cc$trf<-(biomass_all_cc$tn*biomass_all_cc$recollection_female)

boxplot(per_capita_female_biomass~Treatment, data=biomass_all_cc)


#reduced model
biomass_red_nlme<-lme(per_capita_female_biomass~Treatment*Year+recollection_female,
                   random = ~1|trf,biomass_all_cc)
hist(resid(biomass_red_nlme))
anova(biomass_red_nlme)


#plotting ANCOVA figure for RAW data; this is what I present in the manuscript####

#Order the treatments: Low --> Med --> High
repro$Treatment.ord<-ordered(repro$Treatment,levels=c("Low","Medium","High"))

#png("Output/2020.4.11.gobies.per.reef.9.5x5.5.300dpi.png", width = 9.5, height = 5.7, units = 'in', res = 300)

anc1<-ggplot(repro, aes(avg_den_max_seen, egg.week, shape=Treatment.ord, linetype=Treatment.ord, col=Treatment.ord)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x=(bquote('Gobies '~Reef^ -1*'')),y=(bquote('Reproducton '~Reef^ -1~ Week^-1*'')))+
  expand_limits(y=0)+
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=16)) +
  theme(legend.title =element_text(size=17))+
  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,12000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
anc1

#dev.off()

# rerunning same model, but with trial as fixed factor ####

fixed_T<-aov(egg.week ~ Treatment*Year*avg.inhab*(Trial%in%Year), repro)
hist(resid(fixed_T))
qqnorm(resid(fixed_T))
qqline(resid(fixed_T))
plot(fixed_T)
summary(fixed_T)
anova(fixed_T)

boxplot(egg.week~Treatment,repro)

me2<-update(fixed_T, .~. -(Treatment:Year:avg.inhab:Trial))
summary(me2)
anova(me2)

anova(fixed_T,me2)

me3<-update(me3, .~. -(Year:avg.inhab:Trial))
summary(me3)
anova(me3)

anova(me2,me3) #error?

me4<-update(me3, .~. -(Treatment:Year:avg.inhab))
summary(me4)
anova(me4)

anova(me3,me4)

me5<-update(me4, .~. -(Year:avg.inhab))
summary(me5)
anova(me5)

anova(me4,me5)

m1<-aov(egg.week ~ Treatment*Year+avg.inhab+(Trial%in%Year), repro)
hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
plot(m1)
summary(m1)
anova(m1)


m1a<-aov(egg.week ~ Treatment*Year+avg.inhab+(Trial%in%Year)+Treatment*(Trial%in%Year), repro)
hist(resid(m1a))
qqnorm(resid(m1a))
qqline(resid(m1a))
plot(m1a)
summary(m1a)
anova(m1a)

anova(m1,m1a)

boxplot(egg.week~(Trial%in%Year), repro)

boxplot(egg.week~Treatment*(Trial%in%Year), repro)

# runnning per capita output as response variable with same model
head(repro1)

mcap<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female+(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

boxplot(reproduction_per_capita_female_biomass_per_day_per_gram~(Trial%in%Year), repro1)

mcap1<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female+Treatment*(Trial%in%Year), repro1)
hist(resid(mcap1))
qqnorm(resid(mcap1))
qqline(resid(mcap1))
plot(mcap1)
summary(mcap1)
anova(mcap1) #no trial x treatment effect

boxplot(recollection_female~Treatment, repro1)

mcap2<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female
           +Treatment*(Trial%in%Year)+Treatment*recollection_female, repro1)
hist(resid(mcap2))
qqnorm(resid(mcap2))
qqline(resid(mcap2))
plot(mcap2)
summary(mcap2)
anova(mcap2)

# findings: 1) sig. differences by trial when included as a nested fixed factor. 
# 2) Same qualitative results for per capita reproduction per female biomass per day (reproduction per female per gram per day)
# 3) 

#trial as random
mcap3<-lmer(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female+
              Treatment*recollection_female+(1|Trial)+(1|Treatment:Trial)+(1|recollection_female:Trial), REML = F, repro1)
hist(resid(mcap3))
qqnorm(resid(mcap3))
qqline(resid(mcap3))
plot(mcap3)
summary(mcap3)
anova(mcap3)

### biomass recollected

mcap<-aov(per_capita_female_biomass ~ Treatment*Year+recollection_female+(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

mcap<-aov(per_capita_female_biomass ~ Treatment*Year+recollection_female+Treatment*(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

boxplot(per_capita_female_biomass~Treatment, repro1)

mcap<-lmer(per_capita_female_biomass ~ Treatment*Year+recollection_female+(1|Trial), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

# 8.23.20 power analysis #####

#importing dataset where med and high-risk treatments have been combined (risk) vs. low-risk treatment
# "no_risk"

pa<-read.csv("Data/8_10_23_power_analysis.csv") #includes biomass, densities, and egg counts
head(pa)

pa$avg.inhab<-(ceiling((pa$recollection_male_and_female+20)/2))

pa$sqrt.egg.week<-sqrt(pa$egg.week)

# model with replicates at reef level, accounting for effect of covariate, but nothing else

mod_reef<-lm(sqrt.egg.week~Treatment_power*avg.inhab, data = pa)
hist(resid(mod_reef))
anova(mod_reef)

#dropping interaction with covariate
mod_reef1<-lm(sqrt.egg.week~Treatment_power+avg.inhab, data = pa)
hist(resid(mod_reef1))
anova(mod_reef1)
summary(mod_reef1)

# seeing if sqrt egg counts makes the most sense (this looks better than the square-root transformed data...)
mod_reef_reg<-lm(egg.week~Treatment_power+avg.inhab, data = pa)
hist(resid(mod_reef_reg))
anova(mod_reef_reg)
summary(mod_reef_reg)

emmeans(mod_reef1, pairwise~Treatment_power) #plugged these into the calculation for power 

# model with replicates at trial level, i.e. grouped by trial

data_trials<-with(pa, aggregate((sqrt.egg.week), list(Treatment_power=Treatment_power, Trial=Trial), mean))
data_trials

#exporting this table
write.csv(data_trials,"Data\\power_by_trial.csv", row.names = FALSE)

data_trials1<-read.csv("Data/power_by_trial_density_added.csv") #includes average number of fish recollected
data_trials1$avg.inhab<-(ceiling((data_trials1$recollection+20)/2))
head(data_trials1)

data_trials1$sqrt.egg.week<-sqrt(data_trials1$egg.week)

mod_trial<-lm(sqrt.egg.week~Treatment_power+avg.inhab, data = data_trials1)
hist(resid(mod_trial))
anova(mod_trial)

emmeans(mod_trial, pairwise~Treatment_power) #plugged these into the calculation for power

emmeans(mes8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(sqrt.egg.week~Treatment,data=repro)

# want to see plots of egg vs. density, faceted by trial ####

# faceted by Trial*Year

p<-ggplot(repro, aes(x=avg.inhab, y=egg.week, colour=Treatment)) +
  geom_point(size = 1.0, na.rm= TRUE)+
  stat_smooth(method = lm, size = 0.75, na.rm = TRUE, se = F)+
  theme_bw()+
  theme_classic(base_size = 10)+
  labs(x = "density",
       y = "eggs per week",
       colour = "treatment")+
  facet_wrap(. ~ Trial*Year, nrow = 2)+
  scale_colour_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme(legend.position = "right")
p

# faceted only by year

p<-ggplot(repro, aes(x=avg.inhab, y=egg.week, colour=Treatment)) +
  geom_point(size = 1.0, na.rm= TRUE)+
  stat_smooth(method = lm, size = 0.75, na.rm = TRUE, se = T)+
  theme_bw()+
  theme_classic(base_size = 10)+
  labs(x = "density",
       y = "eggs per week",
       colour = "treatment")+
  facet_wrap(. ~ Year, nrow = 2)+
  scale_colour_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme(legend.position = "right")
p

# viz with sqrt.egg count (not super different)

# year*trial

p<-ggplot(repro, aes(x=avg.inhab, y=sqrt(egg.week), colour=Treatment)) +
  geom_point(size = 1.0, na.rm= TRUE)+
  stat_smooth(method = lm, size = 0.75, na.rm = TRUE, se = F)+
  theme_bw()+
  theme_classic(base_size = 10)+
  labs(x = "density",
       y = "eggs per week",
       colour = "treatment")+
  facet_wrap(. ~ Trial*Year, nrow = 2)+
  scale_colour_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme(legend.position = "right")
p

# just year

p<-ggplot(repro, aes(x=avg.inhab, y=sqrt(egg.week), colour=Treatment)) +
  geom_point(size = 1.0, na.rm= TRUE)+
  stat_smooth(method = lm, size = 0.75, na.rm = TRUE, se = T)+
  theme_bw()+
  theme_classic(base_size = 10)+
  labs(x = "density",
       y = "eggs per week",
       colour = "treatment")+
  facet_wrap(. ~ Year, nrow = 2)+
  scale_colour_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme(legend.position = "right")
p

