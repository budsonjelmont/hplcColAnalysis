# Script to determine normalization factors for Zhuo's data based on the median peak area of
# all the NONphosphorylated peptides recovered. Compiles a data frame of all the peptides 
# common to the 3 HPLC configuration heatmaps, then determines which ones were either
# A.) quantified in all 30 LC/MS runs OR B.) had MS/MS triggered in all the runs.
# Then calculates normalization factors for each expt and generates some summary plots.

library(gdata)
library(ggplot2)

setwd('C:\\Users\\Judson\\Documents\\Long column paper\\rebuilt heatmaps all stds set to 1')

sampleList=list.files()

c15b3=read.xls('JE6_1FDR_15cm_3um_3hr_5replicate_AllStdsSetTo1 FwdDb.xls')
c50b3=read.xls('JE6_1FDR_50cm_3um_3hr_5replicate_AllStdsSetTo1 fwdDb.xls')
c50b19=read.xls('JE6_1FDR_50cm_1_9um_3hr_3replicate_AllStdsSetTo1 fwdDb.xls')

#make unique PSM id
c15b3$pepid=paste(c15b3$peptide.data.2..assigned.sequence, c15b3$'charge.state.peptide',sep='_')
c50b3$pepid=paste(c50b3$peptide.data.2..assigned.sequence, c50b3$'charge.state.peptide',sep='_')
c50b19$pepid=paste(c50b19$peptide.data.2..assigned.sequence, c50b19$'charge.state.peptide',sep='_')

#set phosphopeptide logical column to identify rows where PSM represents a phosphorylated peptide
c15b3$isphospho = FALSE
c15b3$isphospho[grep('\\*', c15b3$peptide.data.2..assigned.sequence)] = TRUE
c50b3$isphospho = FALSE
c50b3$isphospho[grep('\\*', c50b3$peptide.data.2..assigned.sequence)] = TRUE
c50b19$isphospho = FALSE
c50b19$isphospho[grep('\\*', c50b19$peptide.data.2..assigned.sequence)] = TRUE

#make data frame of all PSMs observed
pepids=c(c15b3$pepid,c50b3$pepid,c50b19$pepid)
pepids=pepids[-which(duplicated(pepids))]

#function to get relevant PSM data if it was detected and compile it into a data frame
makeDF=function(x,d){
	if(x %in% d$pepid){
		i=which(d$pepid==x)
		counter=d$counter[i]
		pepid=d$pepid[i]
		isphospho=d$isphospho[i]
		ratio=d$SILAC.ratio.23.for.user.selected.SILAC.timepoint1[i]
		l2ratio=log(ratio,base=2)
		qval=d$qvalues.for.SILAC.timepoint1[i]
		unstim.r1=d$peakarea.manual.1.rep1.thresholded.timepoint1[i]
		unstim.r2=d$peakarea.manual.1.rep2.thresholded.timepoint1[i]
		unstim.r3=d$peakarea.manual.1.rep3.thresholded.timepoint1[i]
		unstim.r4=d$peakarea.manual.1.rep4.thresholded.timepoint1[i]
		unstim.r5=d$peakarea.manual.1.rep5.thresholded.timepoint1[i]
		okt34.r1=d$peakarea.manual.1.rep1.thresholded.timepoint2[i]
		okt34.r2=d$peakarea.manual.1.rep2.thresholded.timepoint2[i]
		okt34.r3=d$peakarea.manual.1.rep3.thresholded.timepoint2[i]
		okt34.r4=d$peakarea.manual.1.rep4.thresholded.timepoint2[i]
		okt34.r5=d$peakarea.manual.1.rep5.thresholded.timepoint2[i]
		unstim.nreps=sum(!is.na(c(unstim.r1,unstim.r2,unstim.r3,unstim.r4,unstim.r5)))
		okt34.nreps=sum(!is.na(c(okt34.r1,okt34.r2,okt34.r3,okt34.r4,okt34.r5)))
		nreps=unstim.nreps+okt34.nreps
		nms2=sum(!is.na(c(d$RTs.manual.1.rep.1.timepoint1[i],d$RTs.manual.1.rep.1.timepoint2[i],
		                d$RTs.manual.1.rep.2.timepoint1[i],d$RTs.manual.1.rep.2.timepoint2[i],
		                d$RTs.manual.1.rep.3.timepoint1[i],d$RTs.manual.1.rep.3.timepoint2[i],
		                d$RTs.manual.1.rep.4.timepoint1[i],d$RTs.manual.1.rep.4.timepoint2[i],
		                d$RTs.manual.1.rep.5.timepoint1[i],d$RTs.manual.1.rep.5.timepoint2[i])))
		return(data.frame(counter,pepid,isphospho,ratio,l2ratio,qval,
			unstim.r1,unstim.r2,unstim.r3,unstim.r4,unstim.r5,
			okt34.r1,okt34.r2,okt34.r3,okt34.r4,okt34.r5,
			unstim.nreps,okt34.nreps,nreps,nms2))
	} else {
		return(c(NA,NA,NA))
	}
}

c15b3dat=lapply(pepids,makeDF,c15b3)
c50b3dat=lapply(pepids,makeDF,c50b3)
c50b19dat=lapply(pepids,makeDF,c50b19)

names(c15b3dat)=names(c50b3dat)=names(c50b19dat)=pepids

c15b3dat=do.call(rbind,c15b3dat)
c50b3dat=do.call(rbind,c50b3dat)
c50b19dat=do.call(rbind,c50b19dat)

#make column names unique before cbind
colnames(c15b3dat)=paste(colnames(c15b3dat),'c15b3',sep='.')
colnames(c50b3dat)=paste(colnames(c50b3dat),'c50b3',sep='.')
colnames(c50b19dat)=paste(colnames(c50b19dat),'c50b19',sep='.')

#make composite data frame
dat=data.frame(cbind(c15b3dat,c50b3dat,c50b19dat))

#Get total number of MS2 calls in row
dat$nms2 = apply(dat,1,function(x){
  return(sum(na.omit(as.numeric(c(x['nms2.c15b3'],x['nms2.c50b3'],x['nms2.c50b19'])))))
})

##################################################################################
###Make CV histograms
##################################################################################

#Function to plot histogram of CVs
tpToCondition = function(tp){
  if(tp==1){
    return('unstim')
  } else if(tp==2){
    return('OKT3/4')
  }
}

vblToMethod = function(str){
  if(str=='c15b3'){
    return('15cm, 3um')
  } else if(str=='c50b3'){
    return('50cm, 3um')
  } else if(str=='c50b19'){
    return('50cm, 1.9um')
  }
}

doCVhist = function(df,tp,methodvbl){
  xstr=paste('peakarea.manual.1.CV.percent.timepoint',tp,sep='')
  xlabstr = paste('CVs of ',tpToCondition(tp),' peak areas',sep='')
  gp=ggplot(df, aes_string(xstr,fill = 'isphospho')) +
    geom_histogram(bins = 500) +
    scale_fill_manual(labels = c('Nonphospho', 'Phospho'), values = c('blue', 'red')) +
    ggtitle(paste('Peak area CVs in',vblToMethod(methodvbl),sep=' ')) +
    xlab(xlabstr) +
    guides(fill=guide_legend('')) +
    theme_bw()
  ggsave(paste('Phospho vs unphospho peak area CVs in',vblToMethod(methodvbl),'histogram.png',sep=' '),plot=gp)
}

doCVdensity = function(df,tp,methodvbl){
  xstr=paste('peakarea.manual.1.CV.percent.timepoint',tp,sep='')
  xlabstr = paste('CVs of ',tpToCondition(tp),' peak areas',sep='')
  gp=ggplot(df, aes_string(xstr)) +
    geom_density(df, aes_string(xstr),alpha=0.25) +
    scale_fill_manual(labels = c('Nonphospho', 'Phospho'), values = c('blue', 'red')) +
    ggtitle(paste('Peak area CVs in',vblToMethod(methodvbl),sep=' ')) +
    xlab(xlabstr) +
    guides(fill=guide_legend('')) +
    theme_bw()
  ggsave(paste('Phospho vs unphospho peak area CVs in',vblToMethod(methodvbl),'density plot.png',sep=' '),plot=gp)
}

doCVhist(c15b3,1,'c15b3')
doCVhist(c50b3,1,'c50b3')
doCVhist(c50b19,1,'c50b19')
#doCVdensity(c50b3,1,'c50b3')

##########################################################################################################################
##########################################################################################################################
#Get peptide IDs of nonphospho peptides that are common to all 3 methods (i.e. have peak area in each of the 30 LC/MS runs)
#Alternatively, get only nonphospho peptides where there is an MS/MS sequence in all 30 runs
##########################################################################################################################
##########################################################################################################################

#15cm,3um
c15b3dat_np = c15b3dat[which(c15b3dat$isphospho.c15b3==FALSE),]
#Get all rows where there is quant data in all 10 runs
#c15b3dat_np10r = c15b3dat_np[ which(c15b3dat_np$nreps.c15b3==10), ]
#Get all rows where there is an MS2 call in all 10 runs
c15b3dat_np10r = c15b3dat_np[ which(c15b3dat_np$nms2.c15b3==10), ]

#50cm,3um
c50b3dat_np = c50b3dat[which(c50b3dat$isphospho.c50b3==FALSE),]
#Get all rows where there is quant data for all 10 experiments
#c50b3dat_np10r = c50b3dat_np[ which(c50b3dat_np$nreps.c50b3==10), ]
#Get all rows where there is an MS2 call in all 10 runs
c50b3dat_np10r = c50b3dat_np[ which(c50b3dat_np$nms2.c50b3==10), ]

#50cm,1.9um
c50b19dat_np = c50b19dat[which(c50b19dat$isphospho.c50b19==FALSE),]
#Get all rows where there is quant data for all 10 experiments
#c50b19dat_np10r = c50b19dat_np[ which(c50b19dat_np$nreps.c50b19==10), ]
#Get all rows where there is an MS2 call in all 10 runs
c50b19dat_np10r = c50b19dat_np[ which(c50b19dat_np$nms2.c50b19==10), ]

#Make list of normalization peptides
normalizeIt = intersect(intersect(c15b3dat_np10r$pepid.c15b3,c50b3dat_np10r$pepid.c50b3), c50b19dat_np10r$pepid.c50b19)

Norm15b3 = c15b3dat_np10r[match(normalizeIt, c15b3dat_np10r$pepid.c15b3),]
Norm50b3 = c50b3dat_np10r[match(normalizeIt, c50b3dat_np10r$pepid.c50b3),]
Norm50b19 = c50b19dat_np10r[match(normalizeIt, c50b19dat_np10r$pepid.c50b19),]
#######################################################################################################
#Calculate normalization factors 
#######################################################################################################
###Median of all shared peptides
#15cm,3um
unstim1.r1.153 = median(Norm15b3$unstim.r1.c15b3)
unstim1.r2.153 = median(Norm15b3$unstim.r2.c15b3)
unstim1.r3.153 = median(Norm15b3$unstim.r3.c15b3)
unstim1.r4.153 = median(Norm15b3$unstim.r4.c15b3)
unstim1.r5.153 = median(Norm15b3$unstim.r5.c15b3)
okt34.r1.153 = median(Norm15b3$okt34.r1.c15b3)
okt34.r2.153 = median(Norm15b3$okt34.r2.c15b3)
okt34.r3.153 = median(Norm15b3$okt34.r3.c15b3)
okt34.r4.153 = median(Norm15b3$okt34.r4.c15b3)
okt34.r5.153 = median(Norm15b3$okt34.r5.c15b3)
Normfactors_153 = data.frame(
    area=c(unstim1.r1.153,unstim1.r2.153,unstim1.r3.153,unstim1.r4.153,unstim1.r5.153,
           okt34.r1.153,okt34.r2.153,okt34.r3.153,okt34.r4.153,okt34.r5.153),
    treatment=c(rep('Unstimulated',5),rep('OKT3/4',5))  
  )
#50cm,3um
unstim1.r1.503 = median(Norm50b3$unstim.r1.c50b3)
unstim1.r2.503 = median(Norm50b3$unstim.r2.c50b3)
unstim1.r3.503 = median(Norm50b3$unstim.r3.c50b3)
unstim1.r4.503 = median(Norm50b3$unstim.r4.c50b3)
unstim1.r5.503 = median(Norm50b3$unstim.r5.c50b3)
okt34.r1.503 = median(Norm50b3$okt34.r1.c50b3)
okt34.r2.503 = median(Norm50b3$okt34.r2.c50b3)
okt34.r3.503 = median(Norm50b3$okt34.r3.c50b3)
okt34.r4.503 = median(Norm50b3$okt34.r4.c50b3)
okt34.r5.503 = median(Norm50b3$okt34.r5.c50b3)
Normfactors_503 = data.frame(
  area=c(unstim1.r1.503,unstim1.r2.503,unstim1.r3.503,unstim1.r4.503,unstim1.r5.503,
        okt34.r1.503,okt34.r2.503,okt34.r3.503,okt34.r4.503,okt34.r5.503),
  treatment=c(rep('Unstimulated',5),rep('OKT3/4',5))
)
#50cm,1.9um
unstim1.r1.5019 = median(Norm50b19$unstim.r1.c50b19)
unstim1.r2.5019 = median(Norm50b19$unstim.r2.c50b19)
unstim1.r3.5019 = median(Norm50b19$unstim.r3.c50b19)
unstim1.r4.5019 = median(Norm50b19$unstim.r4.c50b19)
unstim1.r5.5019 = median(Norm50b19$unstim.r5.c50b19)
okt34.r1.5019 = median(Norm50b19$okt34.r1.c50b19)
okt34.r2.5019 = median(Norm50b19$okt34.r2.c50b19)
okt34.r3.5019 = median(Norm50b19$okt34.r3.c50b19)
okt34.r4.5019 = median(Norm50b19$okt34.r4.c50b19)
okt34.r5.5019 = median(Norm50b19$okt34.r5.c50b19)
Normfactors_5019 = data.frame(
  area=c(unstim1.r1.5019,unstim1.r2.5019,unstim1.r3.5019,unstim1.r4.5019,unstim1.r5.5019,
        okt34.r1.5019,okt34.r2.5019,okt34.r3.5019,okt34.r4.5019,okt34.r5.5019),
  treatment=c(rep('Unstimulated',5),rep('OKT3/4',5))
)

#############################################################
#Make nonphospho normalization peptide peak area histograms
#############################################################
makeNormalPepDist = function(df,colnamestr,methodvbl,tp,rep){
  hist(log(as.numeric(df[,eval(colnamestr)]),base=10),breaks=seq(4,12,by=0.1),
       ylim=c(0,15),xlim=c(4,12),
       xlab=paste('Log10 peak area (',vblToMethod(methodvbl),')',sep=''),
       main=paste(tpToCondition(tp),'rep',rep,sep=' '))
}

#15cm,3um
par(mfrow=c(2,5))
makeNormalPepDist(Norm15b3,'unstim.r1.c15b3','c15b3',1,1)
makeNormalPepDist(Norm15b3,'unstim.r2.c15b3','c15b3',1,2)
makeNormalPepDist(Norm15b3,'unstim.r3.c15b3','c15b3',1,3)
makeNormalPepDist(Norm15b3,'unstim.r4.c15b3','c15b3',1,4)
makeNormalPepDist(Norm15b3,'unstim.r5.c15b3','c15b3',1,5)
makeNormalPepDist(Norm15b3,'okt34.r1.c15b3','c15b3',2,1)
makeNormalPepDist(Norm15b3,'okt34.r2.c15b3','c15b3',2,2)
makeNormalPepDist(Norm15b3,'okt34.r3.c15b3','c15b3',2,3)
makeNormalPepDist(Norm15b3,'okt34.r4.c15b3','c15b3',2,4)
makeNormalPepDist(Norm15b3,'okt34.r5.c15b3','c15b3',2,5)

#50cm,3um
par(mfrow=c(2,5))
makeNormalPepDist(Norm50b3,'unstim.r1.c50b3','c50b3',1,1)
makeNormalPepDist(Norm50b3,'unstim.r2.c50b3','c50b3',1,2)
makeNormalPepDist(Norm50b3,'unstim.r3.c50b3','c50b3',1,3)
makeNormalPepDist(Norm50b3,'unstim.r4.c50b3','c50b3',1,4)
makeNormalPepDist(Norm50b3,'unstim.r5.c50b3','c50b3',1,5)
makeNormalPepDist(Norm50b3,'okt34.r1.c50b3','c50b3',2,1)
makeNormalPepDist(Norm50b3,'okt34.r2.c50b3','c50b3',2,2)
makeNormalPepDist(Norm50b3,'okt34.r3.c50b3','c50b3',2,3)
makeNormalPepDist(Norm50b3,'okt34.r4.c50b3','c50b3',2,4)
makeNormalPepDist(Norm50b3,'okt34.r5.c50b3','c50b3',2,5)

#50cm,1.9um
par(mfrow=c(2,5))
makeNormalPepDist(Norm50b19,'unstim.r1.c50b19','c50b19',1,1)
makeNormalPepDist(Norm50b19,'unstim.r2.c50b19','c50b19',1,2)
makeNormalPepDist(Norm50b19,'unstim.r3.c50b19','c50b19',1,3)
makeNormalPepDist(Norm50b19,'unstim.r4.c50b19','c50b19',1,4)
makeNormalPepDist(Norm50b19,'unstim.r5.c50b19','c50b19',1,5)
makeNormalPepDist(Norm50b19,'okt34.r1.c50b19','c50b19',2,1)
makeNormalPepDist(Norm50b19,'okt34.r2.c50b19','c50b19',2,2)
makeNormalPepDist(Norm50b19,'okt34.r3.c50b19','c50b19',2,3)
makeNormalPepDist(Norm50b19,'okt34.r4.c50b19','c50b19',2,4)
makeNormalPepDist(Norm50b19,'okt34.r5.c50b19','c50b19',2,5)

#Boxplot of nonphospho peptide 'stds'
mydat = data.frame(peakarea=c(Normfactors_153$area,Normfactors_503$area,Normfactors_5019$area),
                  treatment=c(as.character(Normfactors_153$treatment),as.character(Normfactors_503$treatment),as.character(Normfactors_5019$treatment)),
                  config=c(rep('15cm,3um',10),rep('50cm,3um',10),rep('50cm,1.9um',10)))

#reorder factors
mydat$config = factor(mydat$config, levels = levels(factor(mydat$config))[c(1,3,2)])
mydat$treatment = factor(mydat$treatment, levels = levels(factor(mydat$treatment))[c(2,1)])

p=ggplot(aes(y = peakarea, x = config, fill = treatment), data = mydat) + 
  geom_boxplot() + 
  scale_fill_brewer(palette="Set1") + 
  ylab("Nonphospho median peak area") +
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

ggsave('Median nonphospho peak area boxplot.png',plot=p)
