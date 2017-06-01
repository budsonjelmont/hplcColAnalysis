# Makes a density plot showing the distributions of peak area CVs among the 5 replicate unstimulated samples
# for all 3 of the HPLC column configurations. Each dataset (column configuration) is plotted as a different
# series.

library(ggplot2)
library(gdata)

setwd("C:/Users/Judson/Documents/Long column paper/Heatmap dumps")

c50b19 = read.xls("JE6_1FDR_phos_50cm_1_9um_3hr_5replicate.xls")
c15b3 = read.xls("JE6_1FDR_phos_15cm_3um_3hr_5replicate_2.xls")
c50b3 = read.xls("JE6_1FDR_phos_50cm_3um_3hr_5replicate.xls")

c50b19$label = "50cm column, 1.9um bead"
c50b19$cv = apply(c50b19,1,function(x){
		s = sd(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE)
		m = mean(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE)
		return(s/m * 100)
	})

c15b3$label = "15cm column, 3.0um bead"
c15b3$cv = apply(c15b3,1,function(x){
		s = sd(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE)
		m = mean(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE)
		return(s/m * 100)
	})

c50b3$label = "50cm column, 3.0um bead"
c50b3$cv = apply(c50b3,1,function(x){
		s = sd(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE)
		m = mean(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE)
		return(s/m * 100)
	})

dat = data.frame(cv=c(c50b19$cv, c15b3$cv, c50b3$cv),
	label=c(c50b19$label, c15b3$label, c50b3$label))

dat = dat[!is.na(dat$cv),]

p <- ggplot(dat, aes(x=cv, fill=label)) + geom_density(alpha=0.2, position="identity") +
#	scale_x_continuous(breaks=seq(0,1,0.2)) +
	labs(x="% coefficient of variation")
ggsave(filename=paste("CV of unstim cells all 3 datasets.jpg",sep=""), plot=p)