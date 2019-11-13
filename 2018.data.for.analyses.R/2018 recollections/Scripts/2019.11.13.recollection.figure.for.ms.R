# Description: figures for MS, black and white, fixed dimesions of plots
#               and made sure that I have right dimesions
# Author: George C Jarvis
# Date: Wed Nov 13 16:05:00 2019
# Notes: going to color columns by the treatment (same as egg count figure),
#       but I don't need a legend
# --------------

#load packages
library(ggplot2)

#read data
reco<-read.csv("Data/2019.10.8.recollection.data.csv", na.strings = "")

reco.means<-with(reco, aggregate((Count), list(Treatment=Treatment), mean))
#reco.means
#now apply the se function to the 4th column [,3]
reco.means$se<-with(reco, aggregate((Count), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
#reco.means
reco.means$Treatment<-ordered(reco.means$Treatment,levels=c("Low","Medium","High"))

#start here when you get back
png("Output/2019.11.13.bw.9.5x5.5.300dpi.png", width = 9.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(reco.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none") + #scale_fill_discrete(limits=c("Low","Medium", "High")) +
  #theme(legend.key.size = unit(1.3,'line')) + 
  #theme(legend.title=element_text(size=34) , legend.text=element_text(size=20)) 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0,0.45))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Gobies Recollected")

dev.off()