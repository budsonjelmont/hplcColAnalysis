# Quick script to read in a file containing pSer and Angio std peak areas
# for all of Zhuo's LC/MS runs, parse the names of the samples to determine 
# their identity, and generate box plots

setwd("C:/Users/Judson/Documents/Long column paper/")

dat=read.table("Zhuo heatmap stds.tab", sep="\t", header=TRUE)

dat$stim = NA
dat$stim[grep("0min",dat$expt)] = "unstimulated"
dat$stim[grep("3min",dat$expt)] = "OKT3/4"
#Reorder factors so that "unstimulated" box comes first
dat$stim = factor(dat$stim, levels = levels(factor(dat$stim))[c(2,1)])

dat$config = NA
dat$config[grep("15cm_3um",dat$expt)] = "15cm,3um"
dat$config[grep("50cm_3um",dat$expt)] = "50cm,3um"
dat$config[grep("50cm_1_9um",dat$expt)] = "50cm,1.9um"
#Reorder factors so that "15cm,3um" plots come first
dat$config = factor(dat$config, levels = levels(factor(dat$config))[c(1,3,2)])

dat$label = paste(dat$config,dat$stim,sep=",")

###Make plot with ggplot2
library(ggplot2)

p=ggplot(aes(y = std1, x = config, fill = stim), data = dat) + 
	geom_boxplot() + 
	scale_fill_brewer(palette="Set1") + 
	ylab("pSer standard peak area") +
	scale_y_continuous(limits=c(1.5E9,6E9),breaks=seq(2E9,6E9,by=1E9),labels=seq(2E9,6E9,by=1E9)) + 
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		#axis.line = element_line(colour = 'black'), #doesn't work in all versions of ggplot2, so set lines individually instead
		axis.line.x = element_line(color='black', size=0.6),
		axis.line.y = element_line(color='black', size=0.6),
		axis.text.x = element_text(colour='black', size=19),
		axis.text.y = element_text(colour='black', size=19),
		axis.title.x = element_blank(),
		axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=26),
		axis.ticks = element_blank(),
		legend.position = c(.85,.94),
		legend.text = element_text(size = 19),
		legend.title=element_blank()
	)
ggsave('pSer std peak area boxplot.png',plot=p)

###Angio std boxplot

p=ggplot(aes(y = std2, x = config, fill = stim), data = dat) + 
	geom_boxplot() + 
	scale_fill_brewer(palette="Set1") + 
	ylab("Angio standard peak area") +
	scale_y_continuous(limits=c(1.5E9,6E9),breaks=seq(2E9,6E9,by=1E9),labels=seq(2E9,6E9,by=1E9)) + 
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		#axis.line = element_line(colour = 'black'), #doesn't work in all versions of ggplot2, so set lines individually instead
		axis.line.x = element_line(color='black', size=0.6),
		axis.line.y = element_line(color='black', size=0.6),
		axis.text.x = element_text(colour='black', size=19),
		axis.text.y = element_text(colour='black', size=19),
		axis.title.x = element_blank(),
		axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=26),
		axis.ticks = element_blank(),
		legend.position = c(.85,.94),
		legend.text = element_text(size = 19),
		legend.title=element_blank()
	)

ggsave('Angio std peak area boxplot.png',plot=p)