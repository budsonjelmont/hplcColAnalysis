library(gdata)
library(qvalue)

c50b19 = read.table("JE6_1FDR_phos_50cm_1_9um_3hr_5replicate/p_q.csv", sep=",",header=TRUE)
c50b19 = c50b19[,!colnames(c50b19) %in% c("X")] 
c50b3 = read.table("JE6_1FDR_phos_50cm_3um_3hr_5replicate/p_q.csv", sep=",",header=TRUE)
c50b3 = c50b3[,!colnames(c50b3) %in% c("X")]
c15b3 = read.table("JE6_1FDR_phos_15cm_3um_3hr_5replicate/p_q.csv", sep=",",header=TRUE)
c15b3 = c15b3[,!colnames(c15b3) %in% c("X")]

#compute qvals with both Benjamini-Hochberg via p.adjust and Storey method with lambda=1
#(which should be equivalent to BH)
c50b19Q = qvalue(c50b19$pvals, lambda=0)
c50b19$bhqvals_qv = c50b19Q$qvalues
c50b19$bhqvals_bh = p.adjust(c50b19$pvals, method="BH")
c50b3Q = qvalue(c50b3$pvals, lambda=0)
c50b3$bhqvals_qv = c50b3Q$qvalues
c50b3$bhqvals_bh = p.adjust(c50b3$pvals, method="BH")
c15b3Q = qvalue(c15b3$pvals, lambda=0)
c15b3$bhqvals_qv = c15b3Q$qvalues
c15b3$bhqvals_bh = p.adjust(c15b3$pvals, method="BH")