#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file, should be sorted according sampleType")
parser$add_argument("-c", "--countDir", type="character", required=TRUE, help="HTSeqGenie output directory")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id")
parser$add_argument("-m", "--method", type="character", default="vst", help="transform method(vst or rlog) [default %(default)s]")
parser$add_argument("-p", "--pair", action="store_true", default=FALSE, help="if specified samples are paired ")
parser$add_argument("-l", "--library", action="store_true", default=FALSE, help="if specified consider libraryType as a confounding factor")
parser$add_argument("-n", "--nCPU", type="integer", default=8,  help="Number of cpu to use [default %(default)s]")
#parser$add_argument("-n", "--add_numbers", action="store_true", default=FALSE, help="Print line number at the beginning of each line [default]")
#parser$add_argument("infile", nargs=1, help="Infile")
args <- parser$parse_args()
designFile <- args$designFile
countDir <- args$countDir
outDir <- args$outDir
sample.id <- args$sample
trans.method <- args$method
#contrastFile <- "contrast.list"

dir.create(outDir,recursive = T,showWarnings = F)

if( file.access(designFile) == -1) {
	stop(sprintf("Specified file ( %s ) does not exist", designFile))
}
if( file.access(countDir) == -1) {
	stop(sprintf("Specified directory ( %s ) does not exist", countDir))
}
if( file.access(outDir) == -1) {
	dir.create(outDir, showWarnings = FALSE)
}

suppressPackageStartupMessages(library("HTSeqGenie"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("ReportingTools"))
suppressPackageStartupMessages(library("hwriter"))
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("Category"))
suppressPackageStartupMessages(library("KEGG.db"))
suppressPackageStartupMessages(library("PFAM.db"))
suppressPackageStartupMessages(library("gage"))
suppressPackageStartupMessages(library("pathview"))

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
## function definition

readCountTable<-function(design,saveDir)
{
    require("HTSeqGenie")
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

## prepare count data and design matrix
myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","character"))
myDesign$sampleType <- factor(myDesign$sampleType,levels=unique(myDesign$sampleType))

###myCountTable<-readCountTable(myDesign,countDir)
if(dir.exists(countDir)) {
    myCountTable <- readCountTable(myDesign,countDir)
    countData  <-  myCountTable[,c(-1,-2)]
}else {
    myCountTable <- read.table(countDir,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    rownames(myCountTable) <- myCountTable[,1]
    countData  <-  myCountTable[,c(-1,-2)]
    countData <- countData[,rownames(myDesign)]
}

#save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
## DESeq2 processing
if(args$pair)
{
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = myDesign, design = ~ patient + sampleType)
}else
{
    if(args$library)
    {
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = myDesign, design = ~ libType + sampleType)
    }else
    {
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = myDesign, design = ~ sampleType)
    }
}

loginfo(sprintf("Setting up multicore (%d) ...",args$nCPU))
register(MulticoreParam(args$nCPU))
dds <- DESeq(dds,parallel=TRUE)
###save.image(file=paste(outDir,"/DESeq2.RData",sep=""))

## QA and visulization
loginfo(sprintf("vst transformation ..."))
if(trans.method == "vst"){
    vsd.blind <- varianceStabilizingTransformation(dds, blind=T)
    vstMat.blind <- assay(vsd.blind)
}else if(trans.method == "rlog"){
    vsd.blind <- rlog(dds)
    vstMat.blind <- assay(vsd.blind)
}else{
    notAllZero <- (rowSums(counts(dds))>0)
    vstMat.blind <- log2(counts(dds,normalized=TRUE)[notAllZero,] + 1)
}
###save.image(file=sprintf("%s/%s.vst.RData",outDir,sample.id))

loginfo(sprintf("QC ..."))
qcAndVizMat(vstMat.blind,myDesign,outDir,colSet=NULL,intgroup=c("sampleType"),ntop=1000,extra=paste0(".",sample.id),sfilter=NULL, gfilter=NULL,complexHeatmap.use=NULL,clonotype.col=NULL,runNMF=FALSE)
	
save.image(file=sprintf("%s/%s.vst.RData",outDir,sample.id))
