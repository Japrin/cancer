#!/usr/bin/env Rscript

args<-commandArgs(T)
if(length(args)<2)
{
	cat("usge depth.plot.R <depth> <cummu depth> <outpng> [xlim:150]\n")
	q()
}
library(RColorBrewer)

afile<-args[1]
bfile<-args[2]
ofile<-args[3]
if(length(args)>3)
{
	xx<-numeric(args[4])
}else
{
	xx<-150
}

a<-read.table(afile,header=T,check.names=F)
b<-read.table(bfile,header=T,check.names=F)

gWidth<-1200
gHeight<-600

J<-length(names(a))-1
if(J==1)
{
	my.col<-"red"
}else if(J==2)
{
	my.col<-c("red","blue")
}else if(J>9)
{
	my.col<-rainbow(J)
	gWidth<-2400
	gHeight<-1200
}else
{
	my.col<-brewer.pal(J,"Set1")
}

png(ofile,gWidth,gHeight)
par(mfrow=c(1,2),mar=c(5,5,4,2))

yy1<-a[,2]*100/sum(a[,2])
yy1_max<-max(yy1)*1.25
legend.txt<-names(a)[2:length(names(a))]
plot(x=a$depth,y=yy1,lwd=3,col=my.col[1],type="p",pch=20,cex.axis=2,cex.lab=2,xlab="sequence depth",ylab="Fraction of bases(%)",xlim=c(0,quantile(a$depth,probs=0.95)),ylim=c(0,yy1_max))
if(J>1)
{
	for(j in 1+(2:J))
	{
		points(x=a$depth,y=a[,j]*100/sum(a[,j]),lwd=3,col=my.col[j-1],pch=20)
		#lines(x=a$depth,y=a[,j]*100/sum(a[,j]),lwd=3,col=my.col[j-1],type="l")
	}
}
legend("topright",legend.txt,bty="n",cex=2,fill=my.col)

J<-length(names(b))-1
legend.txt<-names(a)[2:length(names(a))]
#plot(b$depth,b[,2]*100,type="l",col=my.col[1],lwd=3,cex.axis=2,cex.lab=2,xlab="cumulative sequence depth",ylab="Fraction of bases(%)",xlim=c(0,xx))
plot(b$depth,b[,2]*100,type="l",col=my.col[1],lwd=3,cex.axis=2,cex.lab=2,xlab="cumulative sequence depth",ylab="Fraction of bases(%)",xlim=c(0,quantile(a$depth,probs=0.95)))
if(J>1)
{
	for(j in 1+(2:J))
	{
		lines(b$depth,b[,j]*100,type="l",col=my.col[j-1],lwd=3)
	}
}
legend("topright",legend.txt,bty="n",cex=2,lty=1,lwd=3,col=my.col)

dev.off()
