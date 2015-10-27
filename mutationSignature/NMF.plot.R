#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
if(length(args)<4)
{
	cat("NMF.plot.R <W matrix> <W pdf> <H matrix> <H pdf>\n")
	q()
}

ifile<-args[1]
ofile<-args[2]
ifile2<-args[3]
ofile2<-args[4]

#ifile<-"W.txt"
#ofile<-"W.png"
#ifile2<-"H.txt"
#ofile2<-"H.png"


a<-read.table(ifile,header=T, stringsAsFactors=F, check.names=F)
n<-dim(a)[2]-2

a2<-read.table(ifile2,header=T, stringsAsFactors=F, check.names=F)
n2<-dim(a2)[1]
m2<-dim(a2)[2]-2

library(plotrix)
library(RColorBrewer)

#### plot W matrix ####
WPlot <- function (x,i,n, ...)
{
        ya<-(floor(max(x)/5))*6
		#mycol<-brewer.pal(6,"Dark2")
		mycol<-c(rgb(31,190,240,maxColorValue=255),rgb(35,31,32,maxColorValue=255),rgb(231,39,37,maxColorValue=255),rgb(203,202,202,maxColorValue=255),rgb(161,206,99,maxColorValue=255),rgb(238,200,196,maxColorValue=255))		
		mycol2<-c(rep(mycol[1],16),rep(mycol[2],16),rep(mycol[3],16),rep(mycol[4],16),rep(mycol[5],16),rep(mycol[6],16))
		bb_par<-par(mar=c(2,4,1,2)+0.1)
        at<-barplot(x,col=mycol2,space=c(0.4,0.3),border=NA ,beside=T,ylab="",axisnames=F,main="",cex.axis=1.2,cex.lab=1.2,cex.main=1.2,ylim=c(0,0.2),xpd=NA,...)
        
		abline(h=0,xpd=F,lty=2,col="gray")
		abline(h=0.05,xpd=F,lty=2,col="gray")
		abline(h=0.1,xpd=F,lty=2,col="gray")
		abline(h=0.15,xpd=F,lty=2,col="gray")
		abline(h=0.2,xpd=F,lty=2,col="gray")
		
		txtSig=c("A","B","C","D","E","F","G","H","I","J")
		legend("topright",paste("signature ",txtSig[i],sep=""),cex=1.1,bg=brewer.pal(n,"Set2")[i],text.col="white")
		
		mut<-c("C>A","C>G","C>T","T>A","T>C","T>G")
		if(i==1)
		{
			for(j in c(0:5))
			{
				rect(at[1+16*j],0.29,at[16+16*j],0.31,col=mycol[1+j],border=NA)
				text((at[8+16*j]+at[9+16*j])/2,0.35,labels=mut[j+1],cex=1.2)
			}
		}
		if(i==n)
		{
			#print(a$subtypes)
			staxlab(1,at,a[,1],srt=90,cex=0.6,xpd=NA)
			#staxlab(1,at,a$subtypes,srt=90,cex=0.6,xpd=NA)
		}
        #cat(paste0("i=",i," n=",n,"\n"))
		par(bb_par)
}

pdf(ofile,width=10,height=10)
layout(matrix(c(0,1:n,0)))
a_par<-par(cex.axis=1.0,cex.lab=1.0,xpd=NA,mar=c(6,4,4,2)+0.1)
for(i in 1:n)
{		
	WPlot(a[,i+2]/sum(a[,i+2]),i,n)
}
dev.off()
par(a_par)

#### plot H matrix ####

a2.csum<-apply(a2[,3:(2+m2)],2,sum)
a2.freq<-a2
for(i in 1:n2)
{
	for(j in 3:(2+m2))
	{
		a2.freq[i,j]<-a2.freq[i,j]/a2.csum[j-2]
	}
}

HPlot <- function (x,n, ...)
{
	bb_par<-par(mar=c(6,5,4,2)+0.1)
	#x<-x[,order(x[n,],decreasing=T)]
	x<-x[,order(x[1,],decreasing=T)]
	at<-barplot(x,col=brewer.pal(n,"Set2"),border=FALSE,space=c(0),xlab="",ylab="Contribution of Mutational Signatures",main="",axisnames=F,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,cex.names=1.5)
	staxlab(1,at,colnames(x),srt=45,cex=1.2,xpd=NA)
	txtSig=c("A","B","C","D","E","F","G","H","I","J")
	legend("top",legend=paste("signature",txtSig[c(1:n)],sep=""),bty="n",horiz=T,fill=brewer.pal(n,"Set2"),inset=c(0,-0.15),cex=1.0,xpd=NA)
	par(bb_par)
}

pdf(ofile2,width=8,height=6)
HPlot(as.matrix(a2.freq[,3:(2+m2)]),n2)
dev.off()


#### plot cluser dendrogram based on contribution
#png(paste(ofile2,".dendrogram.png",sep=""),width=800,height=600)
#print(paste(ofile2,".dendrogram.pdf",sep=""))
pdf(paste(ofile2,".dendrogram.pdf",sep=""),width=8,height=6)
plot(hclust(dist(t(a2.freq[,c(-1,-2)]))),xlab="",hang=-1)
dev.off()


