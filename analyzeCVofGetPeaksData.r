library(ggplot2)
library(gdata)

setwd("C:/Users/Judson/Documents/Long column paper/Heatmap dumps")

c50b19 = read.xls("JE6_1FDR_phos_50cm_1_9um_3hr_5replicate.xls")
c15b3 = read.xls("JE6_1FDR_phos_15cm_3um_3hr_5replicate_2.xls")
c50b3 = read.xls("JE6_1FDR_phos_50cm_3um_3hr_5replicate.xls")

###50 cm column, 1.9 um beads

c50b19$GetPeaks.boolean.total.for.Zhuo.data.timepoint1 = as.character(c50b19$GetPeaks.boolean.total.for.Zhuo.data.timepoint1)
c50b19$GetPeaks.boolean.total.for.Zhuo.data.timepoint1[is.na(c50b19$GetPeaks.boolean.total.for.Zhuo.data.timepoint1)] = 0
c50b19$avg = apply(c50b19,1,function(x){
		return(mean(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE))
	})
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

c50b19 = c50b19[!is.na(c50b19$cv),]

p <- ggplot(c50b19, aes(x=cv, fill=GetPeaks.boolean.total.for.Zhuo.data.timepoint1)) + geom_density(alpha=0.2, position="identity") +
#	scale_x_continuous(breaks=seq(0,1,0.2)) +
	labs(x="% coefficient of variation") + 
	guides(fill=guide_legend(title="Total manually defined peaks across 5 replicates"))
ggsave(filename=paste("CV of unstim cells in 50cm1.9um sorted by total # of GetPeaks calls.jpg",sep=""), plot=p)

#Calculate total % of GetPeaks calls for column 1
totalPeaks = sum(!is.na(c50b19$peakarea.manual.1.rep1.thresholded.timepoint1),!is.na(c50b19$peakarea.manual.1.rep2.thresholded.timepoint1),
  !is.na(c50b19$peakarea.manual.1.rep3.thresholded.timepoint1),!is.na(c50b19$peakarea.manual.1.rep4.thresholded.timepoint1),
  !is.na(c50b19$peakarea.manual.1.rep5.thresholded.timepoint1))
totalGetPeaksCalls = sum(as.numeric(c50b19$GetPeaks.boolean.total.for.Zhuo.data.timepoint1))
print(totalGetPeaksCalls/totalPeaks)

###15 cm column, 3.0 um beads

c15b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1 = as.character(c15b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1)
c15b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1[is.na(c15b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1)] = 0
c15b3$avg = apply(c15b3,1,function(x){
		return(mean(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE))
	})
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

c15b3 = c15b3[!is.na(c15b3$cv),]

p <- ggplot(c15b3, aes(x=cv, fill=GetPeaks.boolean.total.for.Zhuo.data.timepoint1)) + geom_density(alpha=0.2, position="identity") +
#	scale_x_continuous(breaks=seq(0,1,0.2)) +
	labs(x="% coefficient of variation") + 
	guides(fill=guide_legend(title="Total manually defined peaks across 5 replicates"))
ggsave(filename=paste("CV of unstim cells in 15cm3.0um sorted by total # of GetPeaks calls.jpg",sep=""), plot=p)

#Calculate total % of GetPeaks calls for column 1
totalPeaks = sum(!is.na(c15b3$peakarea.manual.1.rep1.thresholded.timepoint1),!is.na(c15b3$peakarea.manual.1.rep2.thresholded.timepoint1),
  !is.na(c15b3$peakarea.manual.1.rep3.thresholded.timepoint1),!is.na(c15b3$peakarea.manual.1.rep4.thresholded.timepoint1),
  !is.na(c15b3$peakarea.manual.1.rep5.thresholded.timepoint1))
totalGetPeaksCalls = sum(as.numeric(c15b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1))
print(totalGetPeaksCalls/totalPeaks)

###50 cm column, 3.0 um beads

c50b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1 = as.character(c50b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1)
c50b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1[is.na(c50b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1)] = 0
c50b3$avg = apply(c50b3,1,function(x){
		return(mean(c(as.numeric(x["peakarea.manual.1.rep1.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep2.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep3.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep4.thresholded.timepoint1"]),
			as.numeric(x["peakarea.manual.1.rep5.thresholded.timepoint1"])),na.rm=TRUE))
	})
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

c50b3 = c50b3[!is.na(c50b3$cv),]

p <- ggplot(c50b3, aes(x=cv, fill=GetPeaks.boolean.total.for.Zhuo.data.timepoint1)) + geom_density(alpha=0.2, position="identity") +
#	scale_x_continuous(breaks=seq(0,1,0.2)) +
	labs(x="% coefficient of variation") + 
	guides(fill=guide_legend(title="Total manually defined peaks across 5 replicates"))
ggsave(filename=paste("CV of unstim cells in 50cm3.0um sorted by total # of GetPeaks calls.jpg",sep=""), plot=p)

#Calculate total % of GetPeaks calls for column 1
totalPeaks = sum(!is.na(c50b3$peakarea.manual.1.rep1.thresholded.timepoint1),!is.na(c50b3$peakarea.manual.1.rep2.thresholded.timepoint1),
  !is.na(c50b3$peakarea.manual.1.rep3.thresholded.timepoint1),!is.na(c50b3$peakarea.manual.1.rep4.thresholded.timepoint1),
  !is.na(c50b3$peakarea.manual.1.rep5.thresholded.timepoint1))
totalGetPeaksCalls = sum(as.numeric(c50b3$GetPeaks.boolean.total.for.Zhuo.data.timepoint1))
print(totalGetPeaksCalls/totalPeaks)