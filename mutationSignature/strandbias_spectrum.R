require(plotrix)
mybarplot <- function (x, ...) 
{ 
	ya<-(floor(max(x)/5))*6
	at<-barplot(x,col=c("darkred","skyblue"),space=c(0,0.1),beside=T,ylab="Number of mutations",axisnames=F,main="Strand bias of Mutation",cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(0,ya),...) 
	at<-apply(at,2,mean)
	staxlab(1,at,colnames(x),srt=45,cex=1.3)	
	legend("top",legend=c("Transcribed","Untranscribed"),bty="n",horiz=T,fill=c("darkred","skyblue"),cex=2)
}
setwd("./")
a<-read.table("strandbias_spectrum.txt",header=T,sep="	",row.names=1,check.names=F)
a.matrix=as.matrix(a)
a.matrix
png("Strandbias_spectrum.png",width=800,height=600)
par(mar=c(5,5,4,2)+0.1,xpd=F)
mybarplot(a.matrix)
dev.off()
