# Description: code for analyses and plots for reproduction vs. risk in bluebanded gobies
# Author: George C Jarvis
# Date: Mon Jun 14 12:57:46 2021
# Notes:
# --------------

# loading libraries #####

library(car) # version 3.0-1.0
library(tidyverse) # version 1.3.1
library(lme4) # version 1.1-2.1
#library(lmerTest) # version 3.1-0
#library(ggplot2) # version 3.3.3
library(emmeans) # version 1.6.1
library(sciplot) # version 1.2-0
library(pwr) # version 1.3-0

# import dataset #
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


# analyses ####

#[Egg counts] ####

# 1) analyses for egg counts, pooled across all trials (mixed models)
# - plots for egg counts vs. risk, pooled across trials, adjusted for effect of number of gobies inhabiting the reefs
# - power analysis for egg counts, based on expected 40% difference in reproduction between high- and low-risk environments
# -- in fishes (Gilliam, 2010; Mukherjee et al., 2014)
# 2) t-test for cage effects on reproduction in trial 6 only


# notes

# analyzed mixed models with log-likelihood estimates and chi-square tests.
# started with full model, then removed any NS fixed effects
# Mainly interested in treatment effect, whether reproduction differed among years (proxy for tagging method), whether treatment interacted with year, 
# - how density ('avg.inhab') - along with its interactions with other effects - may have affected reproductive output, while accounting for random variance
# - explained by trial
# remove all NS. interactions with the covariate


## Note: use "REML=F" (maximum likelihood estimates) to compare models with log-likelihood estimates (Pinheiro & Bates, 2000; Bolker et al., 2009)
# model structure for random effects: all slopes are the same (assumed effects are the same among trials), but intercepts are random (magnitude differs)

# set SS to type III

options(contrasts = c("contr.sum","contr.poly"))

# models - egg counts ####

# used sqrt-transformed egg counts per reef per week as response variable

# full model, testing all combinations of random and fixed effects (fit is sigular, trial is the only random effect I can include in the model)

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

# model reduction, after removing the nonsignificant interactions with the random effects

# from here on, reduced higher-order terms at P > 0.05

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

# removing trial, running log-likelihood test for trial effect (df = 1)

mes1a<-lm(sqrt.egg.week ~ Treatment + Year + avg.inhab + (Treatment*Year) + (Treatment*avg.inhab) +
             (Year*avg.inhab) + (Treatment*Year*avg.inhab), repro)
Anova(mes1a, type = "III")
summary(mes1a)

2*(logLik(mes1) - logLik(mes1a)) # Chi2 =  0.67
pchisq(2*(logLik(mes1) - logLik(mes1a)), df = 1, lower.tail=F) # P = 0.41, for df = 1

# keeping trial in the model because it explains some of the residual variance (conditioned on fixed effects)

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
# but that I would caution against sedating gobies in MS222)

# removing treatment*year

mes5<-lmer(sqrt.egg.week ~ Treatment + Year + avg.inhab +
             (1|Trial), REML=F, repro)
Anova(mes5, type = "III")

# results:

# > Anova(mes5, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: sqrt.egg.week
#               Chisq Df Pr(>Chisq)    
# (Intercept)  6.8622  1   0.008804 ** 
# Treatment    2.2600  2   0.323034    
# Year         1.1866  1   0.276018    
# avg.inhab   36.9198  1  1.231e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# remove the year effect (again, worth reporting that there was no statistical difference between the two years = no effect of tagging method on repro)

mes6<-lmer(sqrt.egg.week ~ Treatment + avg.inhab +
             (1|Trial), REML=F, repro)
Anova(mes6, type = "III")

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

# > summary(mes6)
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
# avg.inhab      6.552      1.095 104.186   5.982 3.15e-08 *** # estimate for avg.inhab is here (reported in paper, note: on sqrt-transformed scale)
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#            (Intr) Trtmn1 Trtmn2
# Treatment1 -0.041              
# Treatment2 -0.015 -0.484       
# avg.inhab  -0.984  0.034  0.019

# how do reduced models +/- year compare?

AIC(mes5, mes6)

# result:

# > AIC(mes5, mes6)
# df      AIC
# mes5  7 979.3857
# mes6  6 978.4844

# model without year is slightly better than the model with year in it, but not by much

# egg count LS means with 'emmeans' package ####

# emmeans for the fully-reduced model (mes6)

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

# plot for this:

mes6.emm <- emmeans(mes6, ~ Treatment)
plot(mes6.emm)

# based on LS- means: low > high > med

# how does this compare to plots of the raw data?

bargraph.CI(x.factor = Treatment, response = sqrt.egg.week, main="raw data - sqrt repro/reef/week x risk treatment", 
            xlab="Treatment", ylab="reproduction per reef per week", data = repro)

# looks the same as the plots for emmeans, how about the plots with raw data for egg counts?

bargraph.CI(x.factor = Treatment, response = egg.week, main="raw data - egg.week (not sqrt)/reef/week x risk treatment", 
            xlab="Treatment", ylab="reproduction per reef per week", data = repro)

# raw data, not sqrt transformend, looks different from what LS-means suggest

# decided to go with data (means +/- sem) from emmeans

# barplot for egg counts: back-transformed LS-means #####

# for estimates and std. error for plots: (mean - se, mean, mean + se)^2 (each term +/- sem is squared, because these are square root-transformed means)

# exporting results as a df to make back-transforming means and se easier

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

# exporting figure for reproduction ####

pdf("Output/repro_ms_fig_barplot_adjusted_means_for_density.pdf", width = 8, height = 7)

repro_adjusted_plot

dev.off()

# power analysis for egg counts with 'pwr' package ####

# notes:

# based on the model (sqrt.egg.week ~ Treatment + Year + avg.inhab + (Treatment*Year) + (1|Trial), REML = F, repro)

# - pwr.f2.test is correct test for general linear models
# - using numerator df from model output (Treatment numDF = 2)
# - manipulating denDF to represent conservative (denDF = 8) and liberal (denDF > 8) estimates based on design
# -- the most liberal denDF for my treatment factor is 99, which assumes no difference in treatment effects among trials (most liberal assumption)

# rationale for df (with assistance from Steve Dudgeon):

# For your experiment, I came up with 8 df for the error term for Risk as follows:
#   Since Year is also fixed, the Risk X Year term is fixed, so it is not the error term for Risk, but Risk*Trial(Year) is random, 
#   and should be the error term for Risk. Correct?
#   
#   If so, I calculate for 3 levels for Risk, 3 Trials in each of 2 years a df for this term of
# 2*(3-1)*2 = 8

# usage: pwr.f2.test(u = NULL, v = NULL, f2 = NULL, sig.level = 0.05, power = NULL)

# arguments:

# u =  degrees of freedom for numerator
# v =  degrees of freedom for denominator
# f2 = effect size
# sig.level = Significance level (Type I error probability)
# power = Power of test (1 minus Type II error probability)

#NOTE: the parameter that is left as "NULL" is solved for

pwr.f2.test(u = 2, v = 4, f2 = 0.40, sig.level = 0.05, power = NULL) #17% # mainly for figure

pwr.f2.test(u = 2, v = 8, f2 = 0.40, sig.level = 0.05, power = NULL) #32% # what we did, conservatively, based on our model

pwr.f2.test(u = 2, v = 10, f2 = 0.40, sig.level = 0.05, power = NULL) #40%

pwr.f2.test(u = 2, v = 15, f2 = 0.40, sig.level = 0.05, power = NULL) #57%

pwr.f2.test(u = 2, v = 20, f2 = 0.40, sig.level = 0.05, power = NULL) #71% # if we had done 6 trials per year

pwr.f2.test(u = 2, v = 25, f2 = 0.40, sig.level = 0.05, power = NULL) #81% - would have needed ~ 25 replicates to have reasonable power to detect an effect of 40%

# this amounts to 4 more trials per year, over twice the sampling effort we used in our study

pwr.f2.test(u = 2, v = 30, f2 = 0.40, sig.level = 0.05, power = NULL) #88%

pwr.f2.test(u = 2, v = 35, f2 = 0.40, sig.level = 0.05, power = NULL) #93%

pwr.f2.test(u = 2, v = 40, f2 = 0.40, sig.level = 0.05, power = NULL) #96%

# creating dataframe for basic plot of power

num_df <- c(4,8,10,15,20,25,30,35,40)
pwr <- c(17,32,40,57,71,81,88,93,96)

pwr.df <- data.frame(num_df, pwr)

# plot for power analysis ####

pwr.plot <- ggplot(pwr.df, aes(x=num_df, y=pwr)) + 
  geom_smooth(se = FALSE, size = 1)+
  geom_point(size = 1.5)+
  geom_segment(aes(x = 25, xend=25, y = -Inf, yend=81, linetype = 'Theoretical design (7 trials per year)'), colour = "black", size = 0.8) +
  geom_segment(aes(x= -Inf, xend=25, y = 81, yend = 81,linetype = 'Theoretical design (7 trials per year)'), colour = "black", size = 0.8) +
  geom_segment(aes(x = 8, xend=8, y = -Inf, yend=32, linetype = 'Our design (3 trials per year)'), colour = "black", size = 0.8) +
  geom_segment(aes(x= -Inf, xend=8, y = 32, yend = 32, linetype = 'Our design (3 trials per year)'), colour = "black", size = 0.8) + 
  scale_color_manual("", values = c("red", "black")) + 
  theme_classic() + 
  labs(x="Denominator df", y ="Power") +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  scale_y_continuous(expand = c(0,0), limits = c(10,110), breaks = c(20,40,60,80,100))+
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 14)) +
  theme(
    legend.position = c(0.55, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

# exporting plot  
  
pdf("Output/plot_for_power_analysis.pdf", width = 8, height = 7)

pwr.plot

dev.off()

# t test for trial 6: does reproductive output differ between high- partially caged vs. high- uncaged reefs? ####

t.test(t6_t_test$sqrt_high_cage ,t6_t_test$sqrt_high_uncaged)

# results:

# > t.test(t6_t_test$sqrt_high_cage ,t6_t_test$sqrt_high_uncaged)
# 
# Welch Two Sample t-test
# 
# data:  t6_t_test$sqrt_high_cage and t6_t_test$sqrt_high_uncaged
# t = -0.56836, df = 7.4527, p-value = 0.5865
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -23.42178  14.25430
# sample estimates:
#   mean of x mean of y 
# 35.41706  40.00080 

# t tests suggests no difference in reproduction between the two high-risk treatments

# end repro analyses and viz ####
