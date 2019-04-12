# Description: Kathryn figures for total invert abundance
# Author: George C Jarvis
# Date: Mon Apr 01 22:41:25 2019
# --------------

library(sciplot)
library(ggplot2)
library(extrafont)

invert<-read.csv("Data/Inverts.on.algae.csv")
View(invert)

invert$Species<-ordered(invert$Species, levels=c("SAPA","ZOFA","DIUN","SAHO"))
#sapa, zofa, diun, saho

df<-with(invert, aggregate((Total), list(Species=Species), mean))
df
#now apply the se function to the 4th column [,3]
df$se<-with(invert, aggregate((Total), list(Species=Species), function(x) sd(x)/sqrt(length(x))))[,2]
df

p<- ggplot(df, aes(x=Species, y=x, fill=Species)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  #scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none") + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  theme(legend.title=element_text(size=20) , legend.text=element_text(size=18)) + scale_fill_manual(values=c("#0072B2","#D55E00","#009E73","magenta")) + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(size= 25, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
p + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Species", y="Epifaunal Abundance")
