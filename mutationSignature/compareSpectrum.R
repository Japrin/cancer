#!/usr/local/bin/Rscript

Args <- commandArgs(TRUE)

#for(i in 1:length(args))
#{
#	eval(parse(text=args[[i]]))
#}

if(length(Args)<2)
{
	cat("usage:compareSpectrum.R <spectrum> <outDir>\n")
	q()
}
# infile
# outDir
infile<-Args[1]
outDir<-Args[2]


myDist<-function(x,y)
{
	test.o<-chisq.test(as.numeric(x),p=as.numeric(y),rescale.p=T)
	if(is.na(test.o$p.value))
	{
		#print(x)
		#print(y)
		0
	}else
	{
		if(test.o$p.value<2.2e-16 || test.o$p.value==0)
			pvalue<-2.2e-16  
		else
			pvalue<-test.o$p.value 
		-log10(pvalue)
	}
}
myDistClass<-function(X)
{
	n<-nrow(X)
	matDist<-matrix(as.numeric(rep("na",n*n)),nrow=n)
	for(i in 1:(n-1))
	{
		for(j in 2:n)
		{
			t<-myDist(X[i,],X[j,])
			matDist[i,j]=t
			matDist[j,i]=t
		}	
	}
	as.dist(matDist)
}
afunc<-function(x){ x/sum(x)}
color.map <- function(x) { if (x=="this study") "#FF0000" else "#0000FF" }

library("gplots")

myData<-read.table(infile,header=T,row.names=7,sep="\t",check.names=F)
myData<-t(apply(myData,1,afunc))
donorcolors <- unlist(lapply(row.names(myData), color.map))
#heatmap.2(as.matrix(myData),distfun=myDistClass,Colv=NULL,cexRow=0.1,cexCol=1,density.info="none", trace="none",key=TRUE, symkey=FALSE,col=terrain.colors(100))
setwd(outDir)
png(paste("Spectrum.heatmap.png",sep=""),width=800,height=600)
#par(mar=c(7,6,4,2)+0.1)
heatmap.2(as.matrix(myData),Colv=NULL,cexRow=0.2,cexCol=1.2,density.info="none", trace="none",key=TRUE, symkey=FALSE,col=terrain.colors(100),RowSideColors=donorcolors)
dev.off()




