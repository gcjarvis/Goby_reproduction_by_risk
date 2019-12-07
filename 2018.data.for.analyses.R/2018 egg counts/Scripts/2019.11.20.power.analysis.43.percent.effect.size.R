# Description: calculation for power analysis with 43% reduction in output
# Author: George C Jarvis
# Date: Wed Nov 20 07:46:50 2019
# Notes: this effect size of 43% reduction under high risk environment was taken from (Mukherjee et al. 2014)
# --------------

library(pwr)

pwr.f2.test(u = 2, v = 99.09, f2 = 0.43, sig.level = 0.05, power = NULL)



#parameters:

#u degrees of freedom for numerator
#v degrees of freedomfor denominator
#f2 effect size
#sig.level Significance level (Type I error probability)
#power Power of test (1 minus Type II error probability)


pwr.f2.test(u = NULL, v = NULL, f2 = NULL, sig.level = 0.05, power = NULL)