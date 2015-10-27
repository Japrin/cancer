#!/usr/bin/env Rscript

args<-commandArgs(T)
if(length(args)<4)
{
	cat("run.ABSOLUTE.R <seg data> <snv data> <outDir> <sampleID> [disease]\n")
	q()
}
segFile<-args[1]
snvFile<-args[2]
outDir<-args[3]
sampleID<-args[4]

if(length(args)>4)
{
	diseaseName<-args[5]
}else
{
	diseaseName<-"Esophageal.Cancer"
	#diseaseName<-"Esophageal.squamous"
	#diseaseName<-"BRCA"
}

library(ABSOLUTE)
setwd(outDir)
RunAbsolute(seg.dat.fn=segFile,results.dir=outDir,min.ploidy=0.5,max.ploidy=6,platform="Illumina_WES",max.sigma.h=0.02,copy_num_type="total",sigma.p=0,primary.disease=diseaseName,sample.name=sampleID,max.as.seg.count=1000000,maf.fn=snvFile,min.mut.af=0.05,verbose=T,max.non.clonal=0.5,max.neg.genome=0.5)

load(paste(outDir,"/",sampleID,".ABSOLUTE.RData",sep=""))
write.table(seg.dat$mode.res$mode.tab,file=paste(outDir,"/",sampleID,".ABSOLUTE.solution.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
save.image(file=paste(outDir,"/",sampleID,".ABSOLUTE.run.image.RData",sep=""))
