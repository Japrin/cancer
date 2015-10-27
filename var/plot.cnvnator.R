#!/PROJ/GR/share/Software/R/bin/Rscript

args<-commandArgs(T)
if(length(args)<4)
{
	cat("usage: plot.cnvnator.R <infile> <prefix> <beg> <end> <ylim>\n")
	q()
}
infile<-args[1]
prefix<-args[2]
beg<-args[3]
end<-args[4]

#opt_Ylim<-5
#if(length(args)>4)
#{
#}else
#{
#}
a<-read.table(infile,sep="\t",header=F)
colnames(a)<-c("chrom","beg","end","i","RD_merge","RD_partition","RD_GC1","RD_GC2","RD_GC3","RD_raw")

m<-dim(a)[1]
#f<- (beg<=a$beg & a$end<=end)
rd_raw_summary<-summary(a[,"RD_raw"])
#rd_raw_summary
#png("boxplot.png",width=800,height=600)
#boxplot(a[f,"RD_raw"],col="darkblue")
#dev.off()

CHR<-a$chrom[1]
png(paste(prefix,".png",sep=""),width=800,height=600)

#### version 1
plot(a$beg,a$RD_raw,col="black",type="p",pch=16,cex.lab=1.5,xlab=paste("genome coordinate of ",CHR,sep=""),ylab="RD",ylim=c(0,6*(rd_raw_summary[5]-rd_raw_summary[2])+rd_raw_summary[5]),lwd=1.5,main=prefix)
for(i in 1:m)
{
	lines(c(a[i,"beg"],a[i,"end"]),c(a[i,"RD_merge"],a[i,"RD_merge"]),col="green",lwd=1.5)
}
w=a[1,"end"]-a[1,"beg"]+1
abline(v=as.integer(beg)+w,lwd=1.5,col="darkred")
abline(v=as.integer(end)+w,lwd=1.5,col="darkred")
#print(prefix)
#print(beg)
#print(end)

#### version 2
#plot(a$beg,a$RD_raw,col="black",type="s",pch=16,cex.lab=1.5,xlab=paste("genome coordinate of ",CHR,sep=""),ylab="RD",ylim=c(0,5*(rd_raw_summary[5]-rd_raw_summary[2])+rd_raw_summary[5]),lwd=1.5,main=prefix)
#lines(a$beg,type="s",a$RD_merge,col="green",pch=16,lwd=1.5)
#abline(v=beg,lwd=1.5,col="darkred")
#abline(v=end,lwd=1.5,col="darkred")

dev.off()
