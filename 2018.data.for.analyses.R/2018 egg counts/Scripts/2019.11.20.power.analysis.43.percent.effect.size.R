# Description: calculation for power analysis with 43% reduction in output
# Author: George C Jarvis
# Date: Wed Nov 20 07:46:50 2019
# Notes: this effect size of 43% reduction under high risk environment was taken from (Mukherjee et al. 2014)
# --------------

library(pwr)
library(simr) #testing effect sizes of mixed models

pwr.f2.test(u = 2, v = 99.09, f2 = 0.43, sig.level = 0.05, power = NULL)


#2020.1.8.retesting

pwr.f2.test(u = 2, v = 99, f2 = NULL, sig.level = 0.05, power = .80)

#ran model as linear model, with trial included as a fixed effect

pwr.f2.test(u = 2, v = 86, f2 = NULL, sig.level = 0.05, power = .80)

#still very high power to detect differences based on our sample size (~100% chance)

#parameters:

#u degrees of freedom for numerator
#v degrees of freedom for denominator
#f2 effect size
#sig.level Significance level (Type I error probability)
#power Power of test (1 minus Type II error probability)


pwr.f2.test(u = NULL, v = NULL, f2 = NULL, sig.level = 0.05, power = NULL)

#2020.1.8 update, using fully-reduced model from "2019.12.12.testing.out.l.schuster.code" script

#trying out new package to test for power of mixed models

#SIMR: Power Analysis for Generalised Linear Mixed Models by Simulation
#NOTE: I don't think any of this is right

# I'm not sure how to do a power analysis in R with a mixed model
# I think doing it as a linear model with no random effects is wrong

#1) now that I have the model, I can test the effects, might not work with nlme though...
#idea is that you can change the effect size of the fixed effect

#original model from other script

mod2.2.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")

#ordered treatments so I can see which is which in output
mod2.2.luk<-lme(egg.week~(treatment.ordered*Year.fact)+ treatment.ordered+
                  avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")

#trying to simplify model, just using treatment by trial nested within year
mod2.2.simp<-lme(egg.week~Treatment,random=~1|Trial,repro,method="REML")

levels(repro$treatment.ordered)
levels(repro$Treatment)

fixef(mod2.2.simp)


fixef(mod2.2.simp)["Treatment"]
ranef(mod2.2.luk)
## x
## -0.1148147
fixef(model1)["x"] <- -0.05

fixef(model1)["x"]
## x
## -0.1148147
fixef(model1)["x"] <- -0.05

#power analysis
powerSim(mod2.2.luk)
