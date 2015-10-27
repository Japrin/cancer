#!/usr/bin/env Rscript

library(getopt)
spec=matrix(c(
	"mode",      "m", 2, "character", "sequencing mode; paired (default) or unpaired",
	"parallel",  "p", 2, "integer",   "parallel; default 8",
	"binSize",   "b", 2, "integer",   "binSize; default 30000",
	"bamList",   "i", 1, "character", "bamfile list, required",
	"output",    "o", 1, "character", "output prefix, required",
	"help",      "h", 0, "logical",   "help"
	), ncol=5 , byrow=TRUE )
opt = getopt(spec , opt = commandArgs(TRUE))
if(!is.null(opt$help))
{
	cat(getopt(spec , usage=TRUE))
	q(status=1)
}
if(is.null(opt$bamList) || is.null(opt$output) )
{
	cat(getopt(spec , usage=TRUE))
	q()
}
if(is.null(opt$mode)) { opt$mode<-"paired"; }
if(is.null(opt$parallel)) { opt$parallel<-8; }
if(is.null(opt$binSize)) { opt$binSize<-30000; }

##load(".RData")
bamList<-read.table(opt$bamList,header=F,stringsAsFactors=F)
names(bamList)<-c("sampleID","bam","gender")
head(bamList)
library(cn.mops)
bamDataRanges<-getReadCountsFromBAM(BAMFiles=bamList$bam,sampleNames=bamList$sampleID,refSeqName=c(1:22,"X","Y"),mode=opt$mode,parallel=opt$parallel,WL=opt$binSize)

save.image()
## autosome
f_auto<-(seqnames(bamDataRanges) %in% c(1:22))
X_auto<-bamDataRanges[f_auto,]
res_auto<-cn.mops(X_auto,parallel=8,segAlgorithm="DNAcopy",undo.splits="sdundo",undo.SD=3)
#res_auto<-cn.mops(X_auto,parallel=8,segAlgorithm="DNAcopy")
cnvr<-cnvr(res_auto)
cnvs<-cnvs(res_auto)
if(!(length(cnvr)==0 | length(cnvs)==0))
{
	resCN_auto <- calcIntegerCopyNumbers(res_auto)
}
## X
f_X<-(seqnames(bamDataRanges) %in% c("X"))
X_X<-bamDataRanges[f_X,]
ploidy_X<-ifelse(bamList$gender=="M",1,2)

X_male<-X_X[,ploidy_X==1]
res_Xmale<-haplocn.mops(X_male,parallel=8,segAlgorithm="DNAcopy",undo.splits="sdundo",undo.SD=3)
#res_Xmale<-haplocn.mops(X_male,parallel=8,segAlgorithm="DNAcopy")
cnvr<-cnvr(res_Xmale)
cnvs<-cnvs(res_Xmale)
if(!(length(cnvr)==0 | length(cnvs)==0))
{
	resCN_Xmale<-calcIntegerCopyNumbers(res_Xmale)
}

X_female<-X_X[,ploidy_X==2]
res_Xfemale<-cn.mops(X_female,parallel=8,segAlgorithm="DNAcopy",undo.splits="sdundo",undo.SD=3)
#res_Xfemale<-cn.mops(X_female,parallel=8,segAlgorithm="DNAcopy")
cnvr<-cnvr(res_Xfemale)
cnvs<-cnvs(res_Xfemale)
if(!(length(cnvr)==0 | length(cnvs)==0))
{
	resCN_Xfemale<-calcIntegerCopyNumbers(res_Xfemale)
}
##X_X.norm<-normalizeChromosomes(X_X,ploidy=ploidy_X)
##res_X<-cn.mops(X_X.norm,norm=FALSE,parallel=8,segAlgorithm="DNAcopy",undo.splits="sdundo",undo.SD=3)
##resCN_X <- calcIntegerCopyNumbers(res_X)

## Y

save.image()

write.table(as.data.frame(cnvs(resCN_auto)),file=paste(opt$output,".cnmops.called.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
if(exists("resCN_Xmale"))
{
	write.table(as.data.frame(cnvs(resCN_Xmale)),file=paste(opt$output,".cnmops.called.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F,append=T)
}
if(exists("resCN_Xfemale"))
{
	write.table(as.data.frame(cnvs(resCN_Xfemale)),file=paste(opt$output,".cnmops.called.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F,append=T)
}
save.image()
