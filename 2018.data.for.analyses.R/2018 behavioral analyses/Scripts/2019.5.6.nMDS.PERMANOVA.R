# Description: nMDS for behaviors, with PERMANOVA results
# Author: George C Jarvis
# Date: Mon May 06 11:58:15 2019
# --------------

library(vegan)

behave.2017<-read.csv("Data/2019.5.6..behavior.ss.trial.1.3.csv")
behave.2017<-behave.2017[,3:8]
behave.2018<-read.csv("Data/2019.5.6..behavior.ss.trial.4.5.csv")
behave.2018<-cbind(behave.2018$Trial,behave.2018$Treatment,behave.2018$proportion.exposed,
                   behave.2018$movements.min,behave.2018$bites.min,behave.2018$total.dist.moved,
                   behave.2018$courtship.min)
behave.2018<-behave.2018[,3:8]
behave.2018<-behave.2018[,-c(2,10)]# yay! got rid of the columns I didn't need
View(behave.2018)

#2018 comparison
#create	the	ordination	output	using	bray	curtis
ord<-metaMDS(behave.2018[,-2],k=2,distance='euclidean')
ord$stress

#	if	it	does	not	converge	add	more iterations
#ord<-metaMDS(PercentTotal[,-1],k=2,	distance='bray',	trymax	=	100)	#add	more	iterations
#let's	look	at	the	2D	stress.	Is	it	<	0.3?	
#ord$stress
#	It	is	0.2	which	is	"good/ok".	So,	let's	continue
#	Let's	look	at	the	stress	plot

stressplot(ord)

ordiplot(ord)	

ordiplot(ord,	type	=	'text') #allows you to see which factors are driving trends
View(behave.2018)

#make note of axis limits bacause final plto will have no numbers on the axes
# x axis (-0.2 --> 0.2), y axis (-0.15,0.07)

#	let's	make	a	pretty	plot

png(filename = "Output/nmds.2018.trial.stdev.png", width = 500, height = 300)

plot(1, type='n', xlim=c(-0.15,0.15), ylim=c(-0.075,0.07), xlab='nMDS1', ylab='nMDS2', xaxt='n', yaxt='n')

points(ord$points[behave.2018$Treatment=='Low',1],ord$points[behave.2018$Treatment=='Low',2],	
       pch=20,	col="#0072B2",	cex=2)
points(ord$points[behave.2018$Treatment=='Medium',1],ord$points[behave.2018$Treatment=='Medium',2],	
       pch=20,	col="#009E73",	cex=2)
points(ord$points[behave.2018$Treatment=='High',1],ord$points[behave.2018$Treatment=='High',2],	
       pch=20,	col="#D55E00",	cex=2)

#let's try this by trial
points(ord$points[behave.2018$Trial=='4',1],ord$points[behave.2018$Trial=='4',2],	
       pch=20,	col="#0072B2",	cex=2)
points(ord$points[behave.2018$Trial=='5',1],ord$points[behave.2018$Trial=='5',2],	
       pch=20,	col="#009E73",	cex=2)

#	if	you	want	to	make	the	circles	Standard	deviations
ordiellipse(ord,	groups=behave.2018$Trial,	kind='sd',	border='white',	
            col=c('red','blue'),	lwd=2,	draw	='polygon')

dev.off()

# if you want to outline the outer points
ordihull(ord,	groups=behave.2018$Treatment,	col=c("#009E73","#D55E00","#0072B2"))

#make	a	spider	plot
legend('topleft',	legend	=	paste('2D	stress	=	',	round(ord$stress,2)),	bty='n')
legend('topright',legend=c('Low','Medium','High'),
       col=c("#0072B2","#009E73","#D55E00"),	pch=19,	bty='n')

#if	you	want	to	make	the	circles	Standard	deviations
ordiellipse(ord,	groups=behave.2018$Trial,	kind='sd',	border='white',	
            col=c('palegreen2','purple2'),	lwd=2,	draw	='polygon')

dev.off()

#low, med, high
#"#0072B2","#009E73","#D55E00"

levels(behave.2018.t4.5$Treatment)

#2017 comparison
#create	the	ordination	output	using	bray	curtis
ord<-metaMDS(behave.2017[,-1],k=2,distance='euclidean')

ord$stress #0.11

stressplot(ord)

ordiplot(ord)	

ordiplot(ord,	type	=	'text')

#make note of axis limits bacause final plto will have no numbers on the axes
# x axis (-0.4,0.4), y axis (-0.21,-0.21)

#	let's	make	a	pretty	plot

png(filename = "Output/nmds.2017.png", width = 600, height = 500)

#plot(1, type='n', xlim=c(-0.6,0.6), ylim=c(-0.22,-0.21), xlab='nMDS1', ylab='nMDS2', xaxt='n', yaxt='n')
#plot(1, type='n', xlab='nMDS1', ylab='nMDS2', xaxt='n', yaxt='n')
plot(1, type='n', xlim=c(-0.4,0.4), ylim=c(-0.25,0.25), xlab='nMDS1', ylab='nMDS2', xaxt='n', yaxt='n', cex.axis=1.65)

points(ord$points[behave.2017$Treatment=='Low',1],ord$points[behave.2017$Treatment=='Low',2],	
       pch=20,	col="#0072B2",	cex=2)
points(ord$points[behave.2017$Treatment=='Medium',1],ord$points[behave.2017$Treatment=='Medium',2],	
       pch=20,	col="#009E73",	cex=2)
points(ord$points[behave.2017$Treatment=='High',1],ord$points[behave.2017$Treatment=='High',2],	
       pch=20,	col="#D55E00",	cex=2)

#	if	you	want	to	make	the	circles	Standard	deviations
#ordiellipse(ord,	groups=behave.2018.t4.5$Treatment,	kind='sd',	border='white',	
#            col=c('red','blue','orange'),	lwd=2,	draw	='polygon')

# if you want to outline the outer points by treatment
ordihull(ord,	groups=behave.2018$Treatment,	col=c("#009E73","#D55E00","#0072B2"))
#if you want to outline outer points by trial
ordihull(ord,	groups=behave.2018$Trial,	col=c("#009E73","#D55E00"))

#make	a	spider	plot
legend('topleft',	legend	=	paste('2D	stress	=	',	round(ord$stress,2)),	bty='n')
legend('topright',legend=c('Low','Medium','High'),
       col=c("#0072B2","#009E73","#D55E00"),	pch=19,	bty='n')

dev.off()

#low, med, high
#"#0072B2","#009E73","#D55E00"

levels(behave.2018.t4.5$Treatment)


######making another plot for 2018 behavioral data
png(filename = "Output/nmds.2018.trial.outer.edges.points.png", width = 500, height = 300)

plot(1, type='n', xlim=c(-0.15,0.15), ylim=c(-0.075,0.07), xlab='nMDS1', ylab='nMDS2', xaxt='n', yaxt='n')

points(ord$points[behave.2018$Treatment=='Low',1],ord$points[behave.2018$Treatment=='Low',2],	
       pch=20,	col="#0072B2",	cex=2)
points(ord$points[behave.2018$Treatment=='Medium',1],ord$points[behave.2018$Treatment=='Medium',2],	
       pch=20,	col="#009E73",	cex=2)
points(ord$points[behave.2018$Treatment=='High',1],ord$points[behave.2018$Treatment=='High',2],	
       pch=20,	col="#D55E00",	cex=2)

# if you want to outline the outer points
ordihull(ord,	groups=behave.2018$Treatment,	col=c("#D55E00","#0072B2","#009E73"))

dev.off()