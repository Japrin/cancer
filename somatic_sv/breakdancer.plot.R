#!/usr/local/bin/Rscript

args<-commandArgs(T)
if(length(args)<3)
{
	cat("usage: breakdancer.plot.R <infile> <outDir> <insert size> [read length(default 100)] [n*insert size(default 5)]\n")
	q()
}

#args<-c()
#args[1]<-"test.bed"
#args[2]<-"pic"
#args[3]<-"300"
#args[4]<-"100"
#args[5]<-"5"

library(GenomicFeatures)
library(Gviz)
#library(Homo.sapiens)

infile<-args[1]
outDir<-args[2]
insertSize<-as.numeric(args[3])
if(length(args)>3)
{
	readLen<-as.numeric(args[4])
}else
{
	readLen<-100
}
if(length(args)>4)
{
	nInsertSize<-as.numeric(args[5])
}else
{
	nInsertSize<-5
}

### read in
sv.list<-list()
d<-c()
id<-""
gene<-""
gene_length<-0
sv_size<-0
sv_id<-""
i=0
sv.names<-c()
con <- file(infile,"r")
while(TRUE)
{
	a<-readLines(con,n=1)
	if(is.na(a[1])) { break }
	if(!is.na(grep("^track",a[1],perl=T)[1]))
	{
		if(length(d)>0 && dim(d)[1]>0)
		{
			sv.list[[i]]<-list(data=as.data.frame(d),id=id,gene=gene,gene.length=gene_length,sv.size=sv_size,sv.id=sv_id)
			colnames(sv.list[[i]]$data)<-c("chr","beg","end","name","score","strand","thickStart","thickEnd","itemRgb")
			sv.list[[i]]$data$beg<-as.numeric(gsub("\"","",sv.list[[i]]$data$beg))
			sv.list[[i]]$data$end<-as.numeric(gsub("\"","",sv.list[[i]]$data$end))
			sv.list[[i]]$data$chr<-as.character(gsub("\"","",sv.list[[i]]$data$chr))
			d<-c()
		}
		#gene=SZT2       gene_length=64363       sv_size=343
		m<-regexec("name=(.+?)[[:blank:]]", a[1])
		id<-regmatches(a, m)[[1]][2]
		sv.names<-cbind(sv.names,id)
		m<-regexec("gene=(.+?)[[:blank:]]", a[1])
		gene<-regmatches(a, m)[[1]][2]
		m<-regexec("gene_length=(.+?)[[:blank:]]", a[1])
		gene_length<-as.numeric(regmatches(a, m)[[1]][2])
		m<-regexec("sv_size=([[:digit:]]+)", a[1])
		sv_size<-as.numeric(regmatches(a, m)[[1]][2])
		m<-regexec("SVID=(.+)", a[1])
		sv_id<-as.numeric(regmatches(a, m)[[1]][2])
		i=i+1
	}else
	{
		d<-rbind(d,unlist(strsplit(a[1],"\t",perl=T)))
	}
}
sv.list[[i]]<-list(data=as.data.frame(d),id=id,gene=gene,gene.length=gene_length,sv.size=sv_size,sv.id=sv_id)
colnames(sv.list[[i]]$data)<-c("chr","beg","end","name","score","strand","thickStart","thickEnd","itemRgb")
sv.list[[i]]$data$beg<-as.numeric(gsub("\"","",sv.list[[i]]$data$beg))
sv.list[[i]]$data$end<-as.numeric(gsub("\"","",sv.list[[i]]$data$end))
sv.list[[i]]$data$chr<-as.character(gsub("\"","",sv.list[[i]]$data$chr))
names(sv.list)<-sv.names
close(con)

### refGene prep

refGeneFile <- "/Share/BP/zhenglt/02.pipeline/cancer/database/hg19.refGene.TxDB.sqlite"
refGeneDB <- loadDb(refGeneFile)
###library("TxDb.Hsapiens.UCSC.hg19.knownGene")
###geneModelDB <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
bandFile <- "/Share/BP/zhenglt/02.pipeline/cancer/database/hg19.ideogram.table"
band<-read.table(bandFile,header=T)
### plot ###
plot.sv<-function(sv)
{
	if(sv$gene.length>100000)
	{
		#refGeneFlank<-max(1000*(1.55*as.integer(sv$gene.length/(2*1000))),100000)
		refGeneFlank<-1.3*as.integer(sv$gene.length)
	}else
	{
		refGeneFlank<-100000
	}
	#print(sv$gene.length)
	#print(sv$sv.size)
	#print(refGeneFlank)
	if(!is.na(grep("CTX|ITX",sv$id,perl=T)[1]) | sv$sv.size>100000)
	{
		# translocation
		chr<-sv$data$chr[1:2]
		p1<-(1:length(sv$data$chr))%%2==1
		gstart<-c()
		gend<-c()
		gstart[1]<-min(sv$data$beg[p1])-nInsertSize*insertSize
		gend[1]<-max(sv$data$end[p1])+readLen+nInsertSize*insertSize
		gstart[2]<-min(sv$data$beg[!p1])-nInsertSize*insertSize
		gend[2]<-max(sv$data$end[!p1])+readLen+nInsertSize*insertSize
		
		### ideogram
		ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr[1],band=band, cex=1.5)
		### axis
		axisTrack <- GenomeAxisTrack(cex=1.5)
		### sv track
		svTrack<-list()
		svTrack[[1]]<-AnnotationTrack(start=sv$data$beg[p1],end=sv$data$end[p1],strand=sv$data$strand[p1],chromosome=chr[1],group=sv$data$name[p1],genome="hg19",name="ReadsTrack",shape="arrow",cex=1.5)
		svTrack[[2]]<-AnnotationTrack(start=sv$data$beg[!p1],end=sv$data$end[!p1],strand=sv$data$strand[!p1],chromosome=chr[2],group=sv$data$name[!p1],genome="hg19",name="ReadsTrack",shape="arrow",cex=1.5)
		# detail
		detailTrack<-list()
		detailTrack[[1]] <- GeneRegionTrack(refGeneDB,chromosome=chr[1],genome="hg19",geneSymbols=TRUE,name="refGene",cex.title=1.5,cex=1.5,showExonId=F)
		detailTrack[[2]] <- GeneRegionTrack(refGeneDB,chromosome=chr[2],genome="hg19",geneSymbols=TRUE,name="refGene",cex.title=1.5,cex=1.5,showExonId=F)
		### temp region
		tmpTrack<-list()
		tmpTrack[[1]]<-AnnotationTrack(start=c(gstart[1]),end=c(gend[1]),chromosome=chr[1],genome="hg19",id="tmpId",feature="tmpRegion",group=paste(chr[1],":",gstart[1],"-",gend[1],sep=""),tmpRegion="lightgreen",name="RegionDetail",cex=1.5,details.size=0.75,selectFun=function(identifier, ...) {
			return(TRUE)
		},fun=function(identifier, ...) {
			plotTracks(list(GenomeAxisTrack(scale=0.3,cex=1.5),svTrack[[1]],detailTrack[[1]]),showId=TRUE,add=TRUE,showTitle=FALSE,from=gstart[1],to=gend[1])
		})
		tmpTrack[[2]]<-AnnotationTrack(start=c(gstart[2]),end=c(gend[2]),chromosome=chr[2],genome="hg19",id="tmpId",feature="tmpRegion",group=paste(chr[1],":",gstart[1],"-",gend[1],sep=""),tmpRegion="lightgreen",name="RegionDetail",cex=1.5,details.size=0.75,selectFun=function(identifier, ...) {
			return(TRUE)
		},fun=function(identifier, ...) {
			plotTracks(list(GenomeAxisTrack(scale=0.3,cex=1.5),svTrack[[2]],detailTrack[[2]]),showId=TRUE,add=TRUE,showTitle=FALSE,from=gstart[2],to=gend[2])
		})
		refGeneTrack<-list()
		refGeneTrack[[1]] <- GeneRegionTrack(refGeneDB,chromosome=chr[1],geneSymbols=TRUE,genome="hg19",name="refGene",cex.title=1.5,cex=1.5,showExonId=F)
		refGeneTrack[[2]] <- GeneRegionTrack(refGeneDB,chromosome=chr[2],geneSymbols=TRUE,genome="hg19",name="refGene",cex.title=1.5,cex=1.5,showExonId=F)
		
		### plot
		png(paste(outDir,"/SV_",sv$gene,"_",sv$id,".",sv$sv.id,".png",sep=""),width=1600,height=800)
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(1,2)))

		pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
		plotTracks(list(ideoTrack,axisTrack,tmpTrack[[1]],refGeneTrack[[1]]),groupDetails=TRUE,showId=TRUE,chromosome=chr[1],from=gstart[1]-refGeneFlank,to=gend[1]+refGeneFlank,sizes=c(1,1,8,4), add=TRUE)
		popViewport(1)
		
		pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
		plotTracks(list(ideoTrack,axisTrack,tmpTrack[[2]],refGeneTrack[[2]]),groupDetails=TRUE,showId=TRUE,chromosome=chr[2],from=gstart[2]-refGeneFlank,to=gend[2]+refGeneFlank,sizes=c(1,1,8,4), add=TRUE)
		popViewport(1)

		dev.off()
	}else
	{
		# other
		chr<-sv$data$chr[1]
		gstart<-min(sv$data$beg)-nInsertSize*insertSize
		gend<-max(sv$data$end)+readLen+nInsertSize*insertSize
		### ideogram
		ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr,band=band, cex=1.5)
		### axis
		axisTrack <- GenomeAxisTrack(cex=1.5)
		### sv track
		svTrack<-AnnotationTrack(start=sv$data$beg,end=sv$data$end,strand=sv$data$strand,chromosome=chr,group=sv$data$name,genome="hg19",name="ReadsTrack",from=gstart,to=gend,shape="arrow",cex=1.5)
		### refGene track
		detailRefGeneTrack <- GeneRegionTrack(refGeneDB,chromosome=chr,geneSymbols=TRUE,name="refGene",cex.title=1.5,cex=1.5)
		tmpTrack<-AnnotationTrack(start=c(gstart),end=c(gend),chromosome=chr,genome="hg19",id="tmpId",feature="tmpRegion",group=paste(chr,":",gstart,"-",gend,sep=""),tmpRegion="lightgreen",name="RegionDetail",cex=1.5,details.size=0.75,selectFun=function(identifier, ...) {
			return(identifier=="tmpId")
		},fun=function(identifier, ...) {
			plotTracks(list(GenomeAxisTrack(scale=0.3,cex=1.5),svTrack,detailRefGeneTrack),showId=TRUE,add=TRUE,showTitle=FALSE,from=gstart,to=gend)
		})
		refGeneTrack <- GeneRegionTrack(refGeneDB,chromosome=chr,geneSymbols=TRUE,name="refGene",cex.title=1.5,cex=1.5)
		### plot
		png(paste(outDir,"/SV_",sv$gene,"_",sv$id,".",sv$sv.id,".png",sep=""),width=1600,height=800)
		plotTracks(list(ideoTrack,axisTrack,tmpTrack,refGeneTrack),groupDetails=FALSE,showId=TRUE,from=gstart-refGeneFlank,to=gend+refGeneFlank,sizes=c(1,1,8,4))
		dev.off()
	}
}
ret<-lapply(sv.list,plot.sv)
