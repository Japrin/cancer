#!/usr/local/bin/Rscript

args<-commandArgs(T)
if(length(args)<1)
{
	cat("usage: raw.depth.gatk.R <infile>\n")
	q()
}
infile<-args[1]
#infile<-"/WPS/GR/zhengliangtao/work/BT/BT_21_ganai/pair/stat/cov/RS0130010B/RS0130010B.sample_interval_summary.forPlot.txt"
outDir<-dirname(infile)

#chr     beg     end     total_coverage  average_coverage        RS0130010B_total_cvg    RS0130010B_mean_cvg     RS0130010B_granular_Q1  RS0130010B_granular_median      RS0130010B_granular_Q3  RS0130010B_%_above_1    RS0130010B_%_above_4    RS0130010B_%_above_10   RS0130010B_%_above_15   RS0130010B_%_above_20   RS0130010B_%_above_30   RS0130010B_%_above_40   RS0130010B_%_above_50   RS0130010B_%_above_100

library(ggplot2)

a<-read.table(infile,header=T,stringsAsFactors=F,check.names=F)
a$flag<-rep(".",dim(a)[1])
a$global_bin<-as.numeric(row.names(a))

CHR_LIST<-c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

for(chr in CHR_LIST)
{
	b<-grep(paste("^",chr,"$",sep=""),a$chr,perl=T)
	l<-length(b)
	a[head(b,1),"flag"]<-"tick"
	a[tail(b,1),"flag"]<-"tick"
	a[b[l/2],"flag"]<-"label"
}
f1<-a$flag=="tick"
f2<-a$flag=="label"
if(sum(f2)==0)
{
	f2=f1
}

options(bitmapType='cairo')
png(paste(infile,".raw.depth.gatk.depth.png",sep=""),width=2000,height=600)
plot(a$global_bin,a$average_coverage,type="p",xlab="",ylab="Raw Depth",col="grey",ylim=c(0,100),,xaxt="n",pch=20)
axis(side=1,at=a[f1,"global_bin"],labels=F)
text(a[f2,"global_bin"],par("usr")[3] - 1, srt = 45, adj = 1, labels = a[f2,"chr"], xpd = TRUE)
dev.off()

png(paste(infile,".raw.depth.gatk.coverage.png",sep=""),width=2000,height=600)
plot(a$global_bin,a[,11],type="p",xlab="",ylab="Coverage(%)",col="grey",ylim=c(0,120),,xaxt="n",pch=20)
axis(side=1,at=a[f1,"global_bin"],labels=F)
text(a[f2,"global_bin"],par("usr")[3] - 1, srt = 45, adj = 1, labels = a[f2,"chr"], xpd = TRUE)
dev.off()

png(paste(infile,".raw.depth.gatk.RDdist.png",sep=""),width=800,height=600)
h<-hist(a$average_coverage,breaks=200,plot=F)
plot(h$mids,h$counts,type="p",pch=20,col="darkred",xlab="Average depth per region",ylab="Count",main="Distribution of depth per region")
#plot(h$mids,h$counts/sum(h$counts),type="p",pch=20,col="darkred",xlab="Average depth per region",ylab="",main="Distribution of depth per region")
dev.off()

names(a)<-gsub(".+_above","above",names(a))
a2<-data.frame(y=c(),x=c())
nrow<-dim(a)[1]
ncol<-dim(a)[2]
for(i in 11:(length(a)-2))
{
	a2<-rbind(a2,data.frame(x=rep(names(a)[i],nrow),y=a[,i]))
}
a2$y<-as.numeric(a2$y)

png(paste(infile,".raw.depth.gatk.covdist.png",sep=""),width=800,height=600)
qplot(x,y,data=a2,geom="boxplot",main="Coverage & Depth (per region)")+ geom_point()+xlab("Depth")+ylab("Percentage(%)")
dev.off()
