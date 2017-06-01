# Generates volcano plots and pairwise ratio scatterplots for each of the 3 HPLC column
# configs tested in Zhuo's data.

library(gdata)
library(ggplot2)

setwd('C:\\Users\\Judson\\Documents\\Long column paper\\150522 heatmaps\\after Storey qval recalc')

sampleList=list.files()

c15b3=read.xls('JE6_1FDR_phos_15cm_3um_3hr_5replicate after storey qval recalc.xls')
c50b3=read.xls('JE6_1FDR_phos_50cm_3um_3hr_5replicate after storey qval recalc.xls')
c50b19=read.xls('JE6_1FDR_phos_50cm_1_9um_3hr_5replicate after storey qval recalc.xls')

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

#############################	Pseudo-color scatterplots comparing ratios pairwise between HPLC configs	##############################################
vblToMethod = function(str){
  if(str=='c15b3'){
    return('15cm, 3um')
  } else if(str=='c50b3'){
    return('50cm, 3um')
  } else if(str=='c50b19'){
    return('50cm, 1.9um')
  }
}

makePseudocolorScatter = function(df,xmethodstr,ymethodstr,q01only){
  
  qvalxmethodstr = paste('qval',xmethodstr,sep='.')
  qvalymethodstr = paste('qval',ymethodstr,sep='.')
  
  l2ratioxstr = paste('l2ratio',xmethodstr,sep='.')
  l2ratioystr = paste('l2ratio',ymethodstr,sep='.')
  
  d=df[!is.na(df[,eval(qvalxmethodstr)]),]
  d=d[!is.na(d[,eval(qvalymethodstr)]),]
  
  #filter down by q value
  if(q01only){
    d=d[intersect(which(d[,eval(qvalxmethodstr)]<.01),which(d[,eval(qvalymethodstr)]<.01)),]
  }
  
  #make density color palette
  d$color=densCols(d[,eval(l2ratioystr)],d[,eval(l2ratioxstr)], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  
  #fit line to data and get coefficients to annotate plot
  line = lm(d[,eval(l2ratioystr)]~d[,eval(l2ratioxstr)], data = df)
  m = round(coef(line)[2],3)
  b = round(coef(line)[1],3)
  r2 = round(summary(line)$r.squared,3)
  
  if(as.numeric(b) < 0){
    b = paste('-', abs(b),sep=' ')
  } else {
    b = paste('+',b,sep=' ')
  }
  
  ymxblab = paste('y == ',m,'*x ',b,sep='')
  r2lab = paste('R^2 == ',r2,sep='')
  
  gp=ggplot(d, aes_string(x=l2ratioxstr, y=l2ratioystr, color='color')) +
    scale_color_identity() +
    geom_point(alpha=0.3) +
    ggtitle(paste(vblToMethod(ymethodstr),vblToMethod(xmethodstr),sep=' vs. ')) +
    xlab(paste('log2(OKT34/unstimulated)',vblToMethod(xmethodstr),sep=' ')) +
    ylab(paste('log2(OKT34/unstimulated)',vblToMethod(ymethodstr),sep=' ')) +
    geom_vline(xintercept=0,linetype='dashed') +
    geom_hline(yintercept=0,linetype='dashed') +
    scale_x_continuous(limits=c(-10, 10),labels=seq(-10,10,by=5)) +
    scale_y_continuous(limits=c(-10, 10),labels=seq(-10,10,by=5)) +
    annotate(geom='text', size=c(5,5), x=c(7,7), y=c(9,8), label=c(ymxblab,r2lab), color='blue',parse=TRUE) +
    stat_density2d(aes(fill = ..level..), geom = 'polygon', alpha = 0.25, bins = 150, colour = NA) +
    scale_fill_gradientn(colours=rev(rainbow(100, start=0, end=0.75))) +
    stat_smooth(method = 'lm', formula = y ~ x, col = 'black', linetype='solid', size = 0.5, weight = 2) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(color='black', size=0.6),
      axis.line.y = element_line(color='black', size=0.6),
      axis.text.x = element_text(colour='black', size=12),
      axis.text.y = element_text(colour='black', size=12),
      axis.title.x = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=18),  #ditto
      axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=18),  #ditto
      axis.ticks = element_blank(),
      legend.position = 'none'
    )
  if(q01only){
    ggsave(paste('ratio comparison ',vblToMethod(ymethodstr),' vs ',vblToMethod(xmethodstr),' q01Only pseudocolor.png',sep=''),plot=gp)
  } else {
    ggsave(paste('ratio comparison ',vblToMethod(ymethodstr),' vs ',vblToMethod(xmethodstr),' pseudocolor.png',sep=''),plot=gp)
  }
}

makePseudocolorScatter(dat,'c15b3','c50b3',TRUE)
makePseudocolorScatter(dat,'c15b3','c50b19',TRUE)
makePseudocolorScatter(dat,'c50b3','c50b19',TRUE)

#################################################################################################
#Volcano plots
#################################################################################################
doVolcanoPlot = function(ratios, qvals, methodstr,ratiothresh,qvalthresh){
  
  l2ratios=log(ratios,base=2)
  l10qvals=log(qvals,base=10)
  l2ratiothresh=log(ratiothresh,base=2)
  l10qvalthresh=-1*log(qvalthresh,base=10)
  
  df = data.frame(l2ratios=l2ratios, l10qvals=-1*l10qvals, color=NA)
  df$color[intersect(which(abs(df$l2ratios)>l2ratiothresh),which(df$l2ratios<0))] = 'firebrick2'
  df$color[intersect(which(abs(df$l2ratios)>l2ratiothresh),which(df$l2ratios>0))] = 'green3'
  df$color[which(abs(df$l2ratios)<l2ratiothresh)] = 'black'
  df$color[which(df$l10qvals<l10qvalthresh)] = 'grey'
  
  df$alpha=1
  df$alpha[which(df$l10qvals<l10qvalthresh)]=0
  
  ratiolabelL = paste('<',ratiothresh,' fold',sep='')
  ratiolabelR = paste('>',ratiothresh,' fold',sep='')
  
  gp=ggplot(df, aes(x=l2ratios, y=l10qvals, color=color)) +
    scale_color_identity() +
    geom_point() +
    ggtitle(vblToMethod(methodstr)) +
    xlab('log2(OKT34/unstimulated) peak abundance ratio') +
    ylab('-log10(q value)') +
    geom_vline(xintercept=-1*l2ratiothresh,linetype='dashed',color='blue') +
    geom_vline(xintercept=l2ratiothresh,linetype='dashed',color='blue') +
    geom_hline(yintercept=l10qvalthresh,linetype='dashed',color='blue') +
    scale_x_continuous(limits=c(-10,10),labels=seq(-10,10,by=2),breaks=seq(-10,10,by=2)) +
    scale_y_continuous(limits=c(0,9),labels=seq(0,9,by=1),breaks=seq(0,9,by=1)) +
    annotate(geom='text', size=c(5,5), x=c(-3,3), y=c(9,9), label=c(ratiolabelL,ratiolabelR), color='black',parse=FALSE) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(color='black', size=0.6),
      axis.line.y = element_line(color='black', size=0.6),
      axis.text.x = element_text(colour='black', size=12),
      axis.text.y = element_text(colour='black', size=12),
      axis.title.x = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=18),  #ditto
      axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=18),  #ditto
      axis.ticks = element_blank(),
      legend.position = 'none'
    )

    ggsave(paste(vblToMethod(methodstr),' OG Storeyqvals volcano plot.png',sep=''),plot=gp)
}

doVolcanoPlot(dat$ratio.c15b3,dat$qval.c15b3,'c15b3',2,0.01)
doVolcanoPlot(dat$ratio.c50b3,dat$qval.c50b3,'c50b3',2,0.01)
doVolcanoPlot(dat$ratio.c50b19,dat$qval.c50b19,'c50b19',2,0.01)