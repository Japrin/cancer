#!/usr/bin/Rscript

Args <- commandArgs(TRUE)
infile<-Args[1]
obs_ratio<-as.numeric(Args[2])
fig_file<-Args[3]

if(length(Args)<3)
{
	print("usage:MentoCaloRatoDist.R <infile> <obs_ratio> <fig_file>",quote=F)
	q()
}
out<-c()
a<-read.table(infile)
n=nrow(a)
p<-sum(a[,2]>=obs_ratio)/n
pat=obs_ratio
alpha<-(sort(a[,2],decreasing=T))[0.05*n]

png(fig_file,600,600)
hist(a[,2],breaks=200,col="blue",xlab="NS:S Ratio",main="Distribution of random NS:S ratio")
abline(v=obs_ratio,col="red")
abline(v=alpha,col="red",lty="dashed")
mtext(paste("p-value=",p,sep=""),side=3,at=pat,col="red",line=0)
mtext(sprintf("%4.2f",obs_ratio),side=1,at=obs_ratio,col="red",line=1)
mtext(sprintf("%4.2f",alpha),side=1,at=alpha,col="red",line=1)
