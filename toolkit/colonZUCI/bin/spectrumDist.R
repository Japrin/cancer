#!/usr/bin/Rscript

Args <- commandArgs(TRUE)
infile<-Args[1]
outfile<-Args[2]
iSim<-Args[3]

if(length(Args)<3)
{
	print("usage:spectrumDist.R <infile> <outfile> <simulation times>",quote=F)
	q()
}
out<-c()
a<-read.table(infile)
#G->AnotCpG  G  A GnotCpG.bed.gz  7 8112472
for(i in 1:iSim)
{
	for (j in 1:nrow(a))
	{
		e<-c()
		e<-cbind(e,as.matrix(rep(a[j,4],a[j,5])))	#file
		e<-cbind(e,as.matrix(sample(a[j,6],a[j,5])))	#random line
		e<-cbind(e,as.matrix(rep(a[j,2],a[j,5])))		#ref
		e<-cbind(e,as.matrix(rep(a[j,3],a[j,5])))		#alt
		e<-cbind(e,as.matrix(rep(i,a[j,5])))			#iSim
		e<-cbind(e,as.matrix(rep(a[j,1],a[j,5])))
		e<-cbind(e,as.matrix(rep(a[j,5],a[j,5])))
		e<-cbind(e,as.matrix(rep(a[j,6],a[j,5])))
		out<-rbind(out,e)
	}
}
write.table(out,outfile,sep="\t",row.names=F,col.names=F,quote=F)
