#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-c", "--countDir", type="character", required=TRUE, help="HTSeqGenie output directory")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-t", "--threshold", type="double", default="0.25", help="SIZE FACTOR threshold  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-e", "--ercc", action="store_true", default=FALSE, help="calculate ercc's size factor ? [default %(default)s]")
args <- parser$parse_args()
designFile <- args$designFile
countDir <- args$countDir
out.dir <- args$outDir
sample.id <- args$sample
mode.verbose <- args$verbose
SF.THRESHOLD <- args$threshold
cal.ercc <- args$ercc
print(args)

dir.create(out.dir,recursive=T,showWarnings=F)

#suppressPackageStartupMessages(library(genefilter))
#suppressPackageStartupMessages(library(statmod))
#suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(gplots))
#suppressPackageStartupMessages(library(DESeq2))
#suppressPackageStartupMessages(library(HTSeqGenie))

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
        if(flag==0)
        {
            tEntrez <- as.character(tCount$name)
	    tSymbol <- entrezToXXX(tEntrez)
            out<-data.frame(entrez=tEntrez,symbol=tSymbol)
            #out<-data.frame(entrez=tEntrez)
            flag<-1
        }
        out<-cbind(out,tCount$count)
        
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
    pdf(sprintf("%s.%s",out.prefix,"sizeFactor.dist.pdf"),width = 8,height = 8)
    hist(sf,breaks = 30,freq = F,col = "gray",xlab="Size Factor",main=sprintf("%s",sample.id))
    lines(density(sf),col="orange",lwd=2)
    dev.off()
}

myDesign <- read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))

if(dir.exists(countDir)) {
    myCountTable <- readCountTable(myDesign,countDir)
    countData  <-  myCountTable[,c(-1,-2)]
}else {
    myCountTable <- read.table(countDir,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    rownames(myCountTable) <- myCountTable[,1]
    countData  <-  myCountTable[,c(-1,-2)]
    countData <- countData[,rownames(myDesign)]
}
obj.scdn <- SCDenoise(as.matrix(countData))
obj.scdn <- SCDenoise.normalize(obj.scdn,useERCCSizeFactor = FALSE)

doit <- function(cal.ercc=FALSE)
{
    if(cal.ercc && obj.scdn@withERCC){
        sf <- obj.scdn@size.factor.ERCC
        o.prefix=sprintf("%s/%s.ERCC",out.dir,sample.id)
    }
    else{ 
        sf <- obj.scdn@size.factor.endo
        o.prefix=sprintf("%s/%s.endo",out.dir,sample.id)
    }
    plotSizeFactorDist(sf,o.prefix)
    out.size.factor.df <- data.frame(cellName=names(sf),szieFactor=sf)
    write.table(out.size.factor.df,sprintf("%s.sizeFactor.txt",o.prefix),sep = "\t",row.names = F,col.names = T,quote = F)
    print(quantile(sf,probs=c(0.001,0.002,0.005,0.01,0.05,0.1,0.5,0.9,0.95,0.99,0.995,0.998,0.999)))

    f <- names(sf)[sf > SF.THRESHOLD]
    out.design.df <- data.frame(patient=myDesign$patient,sample=rownames(myDesign))
    out.design.df <- cbind(out.design.df,myDesign[,c(2,3)])
    out.design.df <- out.design.df[f,]
    write.table(out.design.df,sprintf("%s.designSFFiltered.txt",o.prefix),sep = "\t",row.names = F,quote = F)
    return(out.design.df)
}
d.endo.df <- doit(cal.ercc=FALSE)
d.ercc.df <- doit(cal.ercc=TRUE)
f <- intersect(rownames(d.endo.df),rownames(d.ercc.df))
d.df <- d.endo.df[f,,drop=F]
write.table(d.df,sprintf("%s/%s.designSFFiltered.txt",out.dir,sample.id),sep = "\t",row.names = F,quote = F)
save(obj.scdn,file=sprintf("%s/%s.calSizeFactor.RData",out.dir,sample.id))
