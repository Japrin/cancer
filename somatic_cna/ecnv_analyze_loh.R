#!/usr/bin/env Rscript

library(ExomeCNV)

Args <- commandArgs(TRUE)

if(length(Args)<3)
{
	cat("usage:ecnv_analyze_loh.R <normal input> <tumor input> <output>\n")
	q()
}
normal_input<-Args[1]
tumor_input<-Args[2]
output<-Args[3]

setwd(dirname(output))

normal=read.table(normal_input,header=T)
tumor=read.table(tumor_input,header=T)

the.strata=make.loh.strata(normal,tumor)
the.strata.ls=list(the.strata)
#eLOH=LOH.analyze(normal, tumor, strata=the.strata.ls, alpha=0.001, method="deviation.half.norm")
set.seed(50)

### The higher sd.undo is, the coarser the segments will become. The default values are set to 1 and 2.
### The lower alpha is, the coarse the segments will become (harder to accept a breakpoint). The default values are set to 0.05 and 0.01

### optmized parameters for deviation.half.norm:
eLOH=LOH.analyze(normal, tumor, alpha=0.01, method="two.sample.fisher")
the.loh = multi.LOH.analyze(normal, tumor, all.loh.ls=list(eLOH), test.alpha=0.01, method="deviation.half.norm", sdundo=c(3,5), alpha=c(0.01,0.001))

### optmized parameters for variance.f:
### eLOH=LOH.analyze(normal, tumor, alpha=0.05, method="two.sample.fisher")
### the.loh = multi.LOH.analyze(normal, tumor, all.loh.ls=list(eLOH), test.alpha=0.001, method="variance.f", sdundo=c(0,0), alpha=c(0.05,0.01))

#do.plot.loh(the.loh, normal, tumor, "two.sample.fisher", plot.style="baf")
#dev.off()

expanded.loh = expand.loh(the.loh, normal)
write.loh.output(eLOH, paste(output,".eLOH",sep=""))
write.loh.output(the.loh, paste(output,".the",sep=""))
write.loh.output(expanded.loh, paste(output,".expanded",sep=""))

#### admixure estimation
#the.loh[the.loh$LOH==T,]

save.image(paste(output,".RData",sep=""))


