#!/usr/local/bin/Rscript

args<-commandArgs(T)
if(length(args)<3)
{
	cat("usage: cnv-seq.plot.R <infile> <idfile> <outDir> [flank]\n")
	q()
}

library(GenomicFeatures)
library(Gviz)
#library(Homo.sapiens)

### function
getCNV<-function(data, id, upstream = NA, downstream = NA)
{
	sub <- subset(data, cnv == id)
	if (nrow(sub) == 0) 
	{
		return("CNV number given is not avaliable from data")
	}
	if (is.na(upstream)) upstream = 5 * sub$cnv.size[1]
	if (is.na(downstream)) downstream = 5 * sub$cnv.size[1]
	chrom = sub$chromosome[1]
	aa =min(sub$start)
	bb =max(sub$end)
	a = min(sub$start) - upstream
	b = max(sub$end) + downstream
	to.plot <- subset(data, chromosome == chrom & position >= a & position <= b)
	list(cnv.data=to.plot,from=a,to=b,cnv.from=aa,cnv.to=bb)
}


### input

infile<-args[1]
idfile<-args[2]
outDir<-args[3]
if(length(args)>3)
{
	flank=as.numeric(args[4])
}else
{
	flank=as.numeric("NA")
}

data <- read.delim(infile)
idList<-read.table(idfile,header=T,stringsAsFactors=F)

### refGene prep
refGeneFile <- "/WPS/GR/zhengliangtao/00database/UCSC/database/hg19.refGene.withSymbol.sqlite"
refGeneDB <- loadDb(refGeneFile)
bandFile <- "/WPS/GR/zhengliangtao/00database/UCSC/database/hg19.ideogram.table"
band<-read.table(bandFile,header=T)

for(i in 1:dim(idList)[1])
{
	chr <-idList[i,"chr"]
	id  <-idList[i,"id"]
	### ideogram
	ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr,band=band, cex=1.5)
	### axis
	axisTrack <- GenomeAxisTrack(cex=1.5,littleTicks=T)
	### cnv track
	my.cnv.data<-getCNV(data,id,upstream=flank,downstream=flank)
	to.plot<-my.cnv.data$cnv.data
	cnvTrack<-DataTrack(start=to.plot$position, end=to.plot$position, chromosome=chr, genome="hg19", name="log2", data=to.plot$log2,ylim=c(-4,4),cex.axis=1.5,cex.title=1.5,grid=T,col.grid="grey",lty.grid=3,h=10,v=10)

	detailTrack <- GeneRegionTrack(refGeneDB,chromosome=chr,geneSymbols=TRUE,name="refGene",start=my.cnv.data$cnv.from,end=my.cnv.data$cnv.to,cex.title=1.5,cex=1.5)
	tmpTrack<-AnnotationTrack(start=c(my.cnv.data$cnv.from),end=c(my.cnv.data$cnv.to),chromosome=chr,genome="hg19",id="tmpId",feature="tmpRegion",group=paste(chr,":",my.cnv.data$cnv.from,"-",my.cnv.data$cnv.to,sep=""),tmpRegion="lightgreen",name="RegionDetail",details.size=0.7,cex.title=1.5,cex=1.5,fun=function(identifier, ...) {
		plotTracks(list(GenomeAxisTrack(scale=0.3,cex=1.5),detailTrack),showId=TRUE,add=TRUE,showTitle=FALSE)
		})

	### refGene track
	### only specified gene 
	refGeneTrack <- GeneRegionTrack(refGeneDB,chromosome=chr,geneSymbols=TRUE,name="refGene",start=my.cnv.data$from,end=my.cnv.data$to,cex.title=1.5,cex=0.8)
	### plot
	png(paste(outDir,"/",chr,"_",id,".png",sep=""),width=1200,height=800)
	plotTracks(list(ideoTrack,axisTrack,cnvTrack,tmpTrack,refGeneTrack),showId=TRUE,from=my.cnv.data$from,to=my.cnv.data$to,extend.left=1000,extend.right=1000,sizes=c(1,1,10,6,4))
	dev.off()
}
