#!/usr/bin/env Rscript

in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/evaluation/OUT/run01/filtered_TCR_summary/run01.summary.cell.reassigneClonotype.assess.fracAndRecoverRate"
in.file.2 <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/evaluation/liver.P5S4070.tracer.1.3M.N10.txt"
in.file.3 <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/evaluation/OUT/run01/filtered_TCR_summary/run01.summary.cell.reassigneClonotype.assess.tpmAndMinDepth"
out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/evaluation/liver.P5S4070.tracer.1.3M.N10"

in.table <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F)
in.table.2 <-read.delim(in.file.2,header = T,check.names = F,stringsAsFactors = F,sep="\t") 
rownames(in.table.2) <- in.table.2$sample
in.table.3 <-read.delim(in.file.3,header = T,check.names = F,stringsAsFactors = F,sep="\t") 
in.table.3$minReads <- in.table.3$minFrac*in.table.2[in.table.3$CellID,"analyzed"]
in.table.3 <- in.table.3[order(in.table.3$TPM),]

pdf(sprintf("%s.tpmAndMinDepth.pdf",out.prefix),width=8,height = 8)
par(mar=c(5,6,4,4),cex.lab=1.3,cex.axis=1.1,cex.main=1.5)
plot(in.table.3$TPM,in.table.3$minReads,xlab="TPM",ylab="minimum of required reads",log="xy")
abline(h=250000,lty=2,col="red")
dev.off()

#pdf(sprintf("%s.tpmAndMinDepth.example.pdf",out.prefix),width=8,height = 8)
#par(mar=c(5,6,4,4),cex.lab=1.3,cex.axis=1.1,cex.main=1.5)
#d.block <- subset(in.table.3,TPM<50 | minReads>1e6)
#dat.plot <- d.block$minReads
#barplot(dat.plot,xaxt="n")
#axis(1,at = seq_along(dat.plot),labels = d.block$TPM)
#dev.off()

pdf(sprintf("%s.fracAndRecoverRate.pdf",out.prefix),width=8,height = 8)
layout(mat = matrix(c(1:2),nrow=2),heights = c(1:1))
par(mar=c(5,6,4,4),cex.lab=1.3,cex.axis=1.1,cex.main=1.5)
dat.plot.A <- aggregate(isRecoveredA~dFrac,in.table,function(x){ sum(x)/length(x) })
dat.plot.B <- aggregate(isRecoveredB~dFrac,in.table,function(x){ sum(x)/length(x) })
dat.plot <- data.frame(TRA=dat.plot.A$isRecoveredA,TRB=dat.plot.B$isRecoveredB)
rownames(dat.plot) <- sprintf("%d%%",dat.plot.A$dFrac*100)
dat.plot <- t(as.matrix(dat.plot))
barplot(dat.plot,ylim=c(0.75,1),main="Recovery",legend.text=T,args.legend=list(x="topright",inset=c(-0.12,0),xpd=T),xpd=F,beside=T)

dat.plot.nreads <- data.frame()
for(d in unique(in.table$dFrac))
{
    d.block <- data.frame(NReads=in.table.2$analyzed*d,dFrac=d)
    if(nrow(dat.plot.nreads)==0){
        dat.plot.nreads <- d.block
    }else{
        dat.plot.nreads <- rbind(dat.plot.nreads,d.block)
    }
}
#dat.plot.nreads$dFrac <- sprintf("%d%%",dat.plot.nreads$dFrac*100)
#aggregate(NReads~dFrac,dat.plot.nreads,summary)
boxplot(NReads~dFrac,dat.plot.nreads,main="Number of Reads",log="y",xaxt="n")
axis(1,at = seq_along(unique(dat.plot.nreads$dFrac)),labels = sprintf("%d%%",sort(unique(dat.plot.nreads$dFrac)*100)))
dev.off()
