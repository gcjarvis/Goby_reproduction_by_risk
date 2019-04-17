# Description: rethinking df for egg counts
# Author: George C Jarvis
# Date: Tue Apr 16 17:07:25 2019
# --------------

#note: I don't think I should have my data set up by week for my analyses
# if I'm not analyzing by week, I think it inflates my sample size
# and likely artifically inflates my power??

#data wrangling of compiled dataset
egg.den.bio<-read.csv("Data/jarvis.egg.count.data.with.den.max.2019.3.6.csv") #uses adjusted counts for density

#sum of egg counts by treatment, with reef as replicate
#have to do this in excel...