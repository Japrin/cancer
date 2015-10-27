#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)

#for(i in 1:length(args))
#{
#	eval(parse(text=args[[i]]))
#}

if(length(Args)<3)
{
	cat("usage:ecnv_analyze_cnv.R <normal input> <tumor input> <output>\n")
	q()
}

library("ExomeCNV")
library("doMC")

normal_input<-Args[1]
tumor_input<-Args[2]
output<-Args[3]
admixture_rate<-0.3
read_length<-100

setwd(dirname(output))

normal=read.coverage.gatk(normal_input)
tumor=read.coverage.gatk(tumor_input)

e.logR=calculate.logR(normal,tumor)

chr.list = unique(normal$chr)
### exome CNV
registerDoMC()
e.cnv <- foreach(i=1:length(chr.list) , .combine=rbind , .inorder=TRUE) %dopar% {
   idx=which(normal$chr==chr.list[i])
   classify.eCNV(normal=normal[idx,], tumor=tumor[idx,], logR=e.logR[idx], min.spec=0.9999, min.sens=0.9999, option="spec", c=admixture_rate, l=read_length)
}

#e.cnv = classify.eCNV(normal,tumor,logR=e.logR,min.spec=0.9999,min.sens=0.9999,option="spec",c=0.3,l=100)

cnv=multi.CNV.analyze(normal,tumor,logR=e.logR,all.cnv.ls=list(e.cnv),coverage.cutoff=10, min.spec=0.999, min.sens=0.999,option="auc", c=admixture_rate, sdundo=c(1,2), alpha = c(0.05, 0.01) )

do.plot.eCNV(e.cnv,lim.quantile=0.99, style="idx", line.plot=F)
do.plot.eCNV(cnv,lim.quantile=0.99, style="bp",bg.cnv=e.cnv,line.plot=T)
dev.off()

write.output(e.cnv,cnv,paste(output,".cnv",sep=""))

save.image(paste(output,".RData",sep=""))
