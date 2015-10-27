require(plotrix)
library(RColorBrewer)
mybarplot <- function (x, ...) 
{ 
	ya<-(floor(max(x)/5))*6
	at<-barplot(x,col=brewer.pal(6,"Set1"),space=c(0,0.1),beside=T,ylab="Number of mutations",axisnames=F,main="Mutation Spectrum",cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(0,ya),...) 
	at<-apply(at,2,mean)
	staxlab(1,at,colnames(x),srt=45,cex=1.3)
	text(at,x+ya/30,labels=as.character(x))
}
setwd("./")
a<-read.table("spectrum.txt",header=T,sep="\t",check.names=F)
a.matrix=as.matrix(a[,1:6])
a.matrix
png("Spectrum.png",width=800,height=600)
par(mar=c(5,5,4,2)+0.1,xpd=T)
mybarplot(a.matrix)
dev.off()
