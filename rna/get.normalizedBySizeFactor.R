#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-c", "--countDir", type="character", required=TRUE, help="HTSeqGenie output directory")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-e", "--ercc", action="store_true", default=FALSE, help="if specified normalize endo-gene using ERCCs' size factor; otherwise normalize endo-gene using endo-genes' size factor  [default %(default)s]")
parser$add_argument("-n", "--ignoreERCC", action="store_true", default=FALSE, help="whether ignore ERCC data even it's available [default %(default)s]")
parser$add_argument("-x", "--excludeSamples", type="character", help="comma sperated string, specify the samples to be excluded")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
designFile <- args$designFile
countDir <- args$countDir
out.dir <- args$outDir
sample.id <- args$sample
use.ERCC.sf <- args$ercc
mode.verbose <- args$verbose
ignore.ERCC <- args$ignoreERCC
exclude.samples <- args$excludeSamples

print(args)

#### TEST DATA
#designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/marker/P1022.filter.by.marker.design.txt"
#countDir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/OUT/P1022"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/test"
#sample.id <- "P1022"
#do.scLVM <- FALSE
#fit.fdr <- 0.001
#mode.verbose <- FALSE
#cyclebase.rdata <- "/Share/BP/zhenglt/02.pipeline/cancer/rna/data/cycleBase.human.RData"
#use.ERCC.sf <- FALSE
#ignore.ERCC <- TRUE
#sample.id <- "liver"
#designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/sample.design/sample.design.liver.txt"
#countDir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/rc/countDat.liver.txt.gz"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/OUT.ignoreERCC.noTY/liver"


dir.create(out.dir,recursive=T,showWarnings=F)

suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(statmod))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(HTSeqGenie))
suppressPackageStartupMessages(library(scLVM))
suppressPackageStartupMessages(library(rhdf5))

## function definition
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
readCountTable<-function(design,saveDir)
{
    ids<-row.names(design)
    flag<-0
    for(i in 1:length(ids) )
    {
        #tCount <- getTabDataFromFile(paste(countDir,"/",ids[i],sep=""),"counts_gene")
        tCount <- getTabDataFromFile(paste(countDir,"/",ids[i],sep=""),"counts_gene_exonic")
        #print(head(tCount))
        tCount <- tCount[!grepl("^ERCC",tCount$name,perl = T),]
        if(flag==0)
        {
            tEntrez <- as.character(tCount$name)
	        tSymbol <- entrezToXXX(tEntrez)
            out<-data.frame(entrez=tEntrez,symbol=tSymbol)
            #out<-data.frame(entrez=tEntrez)
            flag<-1
        }
        out<-cbind(out,tCount$count)
        cat(sprintf("sample %s done.\n",ids[i]))
        
    }
    names(out)<-c("entrez","symbol",ids)
    rownames(out)<-out$entrez
    return(out)
}

# plot size factor distribution
plotSizeFactorDist <- function(sf,out.prefix,sample.id.toHighlight=NULL)
{
    pdf(sprintf("%s.%s",out.prefix,"sizeFactor.pdf"),width = 10,height = 6)
    par(cex.lab=1.5,mar=c(5,5,4,2)+0.1)
    cdata <- sort(sf)
    ccol <- rep("darkblue",length(sf))
    names(ccol) <- names(cdata)
    if(!is.null(sample.id.toHighlight)) {
        ccol[sample.id.toHighlight] <- "red"
    }
    barplot(cdata,col=ccol,border=NA,xlab="cell index",ylab="size factor",xaxt="n")
    box(which = "inner",lwd=4)
    par(new=TRUE, oma=c(1,7,1,1),mar=c(5,5,1,2)+0.1)
    layout(matrix(4:1,2))
    b.midpoint <- barplot(cdata[1:10],col=ccol[1:10],border=NA,xlab="",ylab="",xaxt="n")
    text(b.midpoint,y=-0.01, srt = 45, adj = 1, labels = names(cdata)[1:10], xpd = TRUE,cex=1.0)
    box(which = "figure")
    dev.off()
}
myDesign <- read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))

if(!is.null(exclude.samples) && exclude.samples != ""){
    if(file.exists(exclude.samples)){
        sample.id.toExclude <- c()
        tryCatch({ 
            sample.id.toExclude <- read.table(exclude.samples,sep = "\t",header = F,check.names = F,stringsAsFactors = F)$V1 
        },error=function(e){e})
    }else{
        sample.id.toExclude <- unlist(strsplit(exclude.samples,split=",",perl=T))
    }
    myDesign <- myDesign[!rownames(myDesign) %in% sample.id.toExclude,]
}

###if(mode.verbose)
###{
out.design.df <- data.frame(patient=myDesign$patient,sample=rownames(myDesign))
out.design.df <- cbind(out.design.df,myDesign[,c(2,3)])
write.table(out.design.df,sprintf("%s/%s.designUsed.txt",out.dir,sample.id),sep = "\t",row.names = F,quote = F)
###}
if(dir.exists(countDir)) {
    myCountTable <- readCountTable(myDesign,countDir)
    countData  <-  myCountTable[,c(-1,-2)]
}else {
    myCountTable <- read.table(countDir,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    rownames(myCountTable) <- myCountTable[,1]
    countData  <-  myCountTable[,c(-1,-2)]
    countData <- countData[,rownames(myDesign)]
}
obj.scdn <- SCDenoise(as.matrix(countData),ignore.ERCC=ignore.ERCC)
obj.scdn <- SCDenoise.normalize(obj.scdn,useERCCSizeFactor = use.ERCC.sf)

#plotSizeFactorDist(obj.scdn@size.factor.endo,sprintf("%s/%s",out.dir,sample.id))

out.df <- data.frame(geneID=rownames(obj.scdn@normalized.endo))
out.df$geneName <- entrezToXXX(out.df$geneID)
out.df <- cbind(out.df, obj.scdn@normalized.endo)
write.table(out.df,file=sprintf("%s/%s.all.countData.sfNormalized.txt",out.dir,sample.id),sep="\t",quote = F,row.names = F)


