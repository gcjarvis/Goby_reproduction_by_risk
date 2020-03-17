# Description: power analysis with 'pwr' package
# Author: George C Jarvis
# Date: Wed Mar 18 09:58:09 2020
# Notes:
# --------------

library(pwr)

# pwr package for power analysis: 
# - pwr.f2.test is correct test for general linear models
# - using numerator df from model output (Treatment numDF = 2)
# - manipulating denDF to represent conservative (denDF = 8) and liberal (denDF > 8) estimates based on model output
# -- the most liberal denDF for my treatment factor is 99, which we all think is incorrect

#parameters:

# u =  degrees of freedom for numerator
# v =  degrees of freedom for denominator
# f2 = effect size
# sig.level = alpha (probability of Type I error)
# power = 1 - beta (1 minus probability of Type II error)

# solving for power (i.e. power parameter is "NULL")

pwr.f2.test(u = 2, v = 8, f2 = 0.40, sig.level = 0.05, power = NULL) #32%

pwr.f2.test(u = 2, v = 10, f2 = 0.40, sig.level = 0.05, power = NULL) #40%

pwr.f2.test(u = 2, v = 20, f2 = 0.40, sig.level = 0.05, power = NULL) #71%

pwr.f2.test(u = 2, v = 25, f2 = 0.40, sig.level = 0.05, power = NULL) #81% - would have needed ~ 25 replicates to detect an effect of 40%

pwr.f2.test(u = 2, v = 30, f2 = 0.40, sig.level = 0.05, power = NULL) #88%

pwr.f2.test(u = 2, v = 40, f2 = 0.40, sig.level = 0.05, power = NULL) #96%

# same analysis, but with effect size found in my paper (~20% between low-risk and medium- + high-risk treatments)

pwr.f2.test(u = 2, v = 8, f2 = 0.20, sig.level = 0.05, power = NULL) #18%

pwr.f2.test(u = 2, v = 10, f2 = 0.20, sig.level = 0.05, power = NULL) #22%

pwr.f2.test(u = 2, v = 20, f2 = 0.20, sig.level = 0.05, power = NULL) #41%

pwr.f2.test(u = 2, v = 30, f2 = 0.20, sig.level = 0.05, power = NULL) #58%

pwr.f2.test(u = 2, v = 40, f2 = 0.20, sig.level = 0.05, power = NULL) #72%

pwr.f2.test(u = 2, v = 50, f2 = 0.20, sig.level = 0.05, power = NULL) #81% #would have needed ~ 50 replicates to detect an effect of 20%

pwr.f2.test(u = 2, v = 60, f2 = 0.20, sig.level = 0.05, power = NULL) #88%

pwr.f2.test(u = 2, v = 70, f2 = 0.20, sig.level = 0.05, power = NULL) #93%

pwr.f2.test(u = 2, v = 80, f2 = 0.20, sig.level = 0.05, power = NULL) #96% #starts to plateau at 80 replicates
