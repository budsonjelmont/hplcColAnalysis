library(ggplot2)
library(gdata)

setwd("C:/Users/Judson/Documents/Long column paper/Heatmap dumps")

c50b19 = read.xls("JE6_1FDR_phos_50cm_1_9um_3hr_5replicate.xls")
c15b3 = read.xls("JE6_1FDR_phos_15cm_3um_3hr_5replicate_2.xls")
c50b3 = read.xls("JE6_1FDR_phos_50cm_3um_3hr_5replicate.xls")

dat_c50b19 = data.frame(peakwidth = rep(NA,nrow(c50b19)*5))
dat_c50b19$peakwidth =c(c50b19$SIC.RT.window.size.manual.1.rep.1.timepoint1,
			c50b19$SIC.RT.window.size.manual.1.rep.2.timepoint1,
			c50b19$SIC.RT.window.size.manual.1.rep.3.timepoint1,
			c50b19$SIC.RT.window.size.manual.1.rep.4.timepoint1,
			c50b19$SIC.RT.window.size.manual.1.rep.5.timepoint1)
dat_c50b19$label = "50cm column, 1.9um bead"

dat_c15b3 = data.frame(peakwidth = rep(NA,nrow(c15b3)*5))
dat_c15b3$peakwidth =c(c15b3$SIC.RT.window.size.manual.1.rep.1.timepoint1,
			c15b3$SIC.RT.window.size.manual.1.rep.2.timepoint1,
			c15b3$SIC.RT.window.size.manual.1.rep.3.timepoint1,
			c15b3$SIC.RT.window.size.manual.1.rep.4.timepoint1,
			c15b3$SIC.RT.window.size.manual.1.rep.5.timepoint1)
dat_c15b3$label = "15cm column, 3.0um bead"

dat_c50b3 = data.frame(peakwidth = rep(NA,nrow(c50b3)*5))
dat_c50b3$peakwidth =c(c50b3$SIC.RT.window.size.manual.1.rep.1.timepoint1,
			c50b3$SIC.RT.window.size.manual.1.rep.2.timepoint1,
			c50b3$SIC.RT.window.size.manual.1.rep.3.timepoint1,
			c50b3$SIC.RT.window.size.manual.1.rep.4.timepoint1,
			c50b3$SIC.RT.window.size.manual.1.rep.5.timepoint1)
dat_c50b3$label = "50cm column, 3.0um bead"

dat = rbind(dat_c50b19, dat_c15b3, dat_c50b3)
dat = dat[!is.na(dat$peakwidth),]

p <- ggplot(dat, aes(x=peakwidth, fill=label)) + geom_density(alpha=0.2, position="identity") +
#	scale_x_continuous(breaks=seq(0,1,0.2)) +
	labs(x="Peak width (s)")
ggsave(filename=paste("Peak width of unstim cells all 3 datasets.jpg",sep=""), plot=p)