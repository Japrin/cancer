#!/usr/bin/env Rscript

args<-commandArgs(T)
if(length(args)<2)
{
	cat("usage: plot.coverage.byChr.R <infile1> <infile2> <outfile>\n")
	q()
}
library(RColorBrewer)

ifile<-args[1]
ifile2<-args[2]
ofile<-args[3]
a<-read.table(ifile,header=T,sep="\t",row.names=1,check.names=F)
a.matrix=as.matrix(a)
b<-read.table(ifile2,header=T,sep="\t",row.names=1,check.names=F)
b.matrix=as.matrix(b)
d<-dim(a.matrix)

gWidth<-1200
gHeight<-600

J<-d[1]
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
#print(a)
#print(b)
png(ofile,gWidth,gHeight)
par(mar=c(5,5,6,6)+0.1,xpd=T)
yMax<-50*ceiling(max(a.matrix)/50)
at<-barplot(a.matrix,col=my.col,beside=T,border=FALSE,xlab="",ylab="Mean depth",main="",ylim=c(0,yMax),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,cex.names=1.5)
legend("top",horiz=TRUE,legend=row.names(a.matrix),bty="n",fill=my.col,inset=c(0,-0.05),cex=1.2)

my.at<-colMeans(at)
par(new=T)
par(mar=c(5,5,6,6)+0.1,xpd=T)
plot(my.at,b[1,],type="b",axes = FALSE, bty = "n", xlab = "", ylab = "",col=my.col[1],ylim=c(0,1))
if(J>1)
{
	for(i in 2:d[1])
	{
		lines(my.at,b[i,],type="b",col=my.col[i])
	}
}
axis(4,ylim=c(0,1),las=1,cex.axis=1.5,cex.lab=1.5)
#axis(4,at=seq(from=0,to=1,by=0.2),ylim=c(0,1),las=1,cex.axis=1.5,cex.lab=1.5)
mtext("Proportion of covered bases",side=4,col="black",line=4,cex=1.5) 
dev.off()
