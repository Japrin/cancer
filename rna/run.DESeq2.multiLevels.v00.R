#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-b", "--contrastFile", type="character", required=TRUE, help="contrast file")
parser$add_argument("-c", "--countDir", type="character", required=TRUE, help="HTSeqGenie output directory")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id")
parser$add_argument("-p", "--pair", action="store_true", default=FALSE, help="if specified samples are paired ")
parser$add_argument("-l", "--library", action="store_true", default=FALSE, help="if specified consider libraryType as a confounding factor")
parser$add_argument("-q", "--qval", default=0.01, type="double", help="qval for differential genes [default %(default)s]")
parser$add_argument("-n", "--nCPU", type="integer", default=8,  help="Number of cpu to use [default %(default)s]")
#parser$add_argument("-n", "--add_numbers", action="store_true", default=FALSE, help="Print line number at the beginning of each line [default]")
#parser$add_argument("infile", nargs=1, help="Infile")
args <- parser$parse_args()
designFile <- args$designFile
countDir <- args$countDir
outDir <- args$outDir
contrastFile <- args$contrastFile
opt.qval <- args$qval
sample.id <- args$sample
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

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")

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

plotCluster<-function(resSigStrict,vsd,designM,outDir)
{
    if(length(levels(designM$sampleType))<=9)
    {
	    require("RColorBrewer")
	    colSet <- brewer.pal(9,"Set1")
	    patientcolors <- colSet[as.numeric(designM$sampleType)]
    }else
    {
	    patientcolors <- as.numeric(designM$sampleType)
    }
    
    pdf(paste(outDir,"/DESeq2.plot.DECluster.pdf",sep=""))
    
    select <- rownames(resSigStrict)[1:30]
    select.f <- !is.na(select)
    select <- select[select.f]
    dat.plot <-assay(vsd)[select,]
    rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
    if(length(select)>0)
    {
    	heatmap.2(dat.plot, ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", trace="none", margin=c(10, 6), main="Most Differentially Expressed genes")
    }
    
    select <- rownames(resSigStrict)
    select.f <- !is.na(select)
    select <- select[select.f]
    dat.plot <-assay(vsd)[select,]
    rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
    if(length(select)>0)
    {
    	heatmap.2(dat.plot, ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", trace="none", margin=c(10, 6), main="All Differentially Expressed genes",labRow=F)
    } 
    dev.off()
}

makePersampleFC<-function(cMat,myDesign)
{
    ##cMat<-vstMat
    fc.df=data.frame(entrez=rownames(cMat),stringsAsFactors=F)
    fc.df$symbol=entrezToXXX(fc.df$entrez)
    fc.df$ensg=entrezToXXX(fc.df$entrez,"ENSG")
    #tumorInfo<-myDesign[myDesign$sampleType=="tumor",]
    #normalInfo<-myDesign[myDesign$sampleType=="normal",]
    tumorInfo <- subset(myDesign,sampleType=="tumor"|sampleType=="treat")
    normalInfo <- subset(myDesign,sampleType=="normal"|sampleType=="control")
    for( i in 1:nrow(tumorInfo))
    {
	    sampleT=rownames(tumorInfo)[i]
	    sampleN=rownames(normalInfo[normalInfo$patient==tumorInfo$patient[i],])
	    sampleFC=cMat[,sampleT]-cMat[,sampleN]
	    fc.df=cbind(fc.df,sampleFC)
    }
    colnames(fc.df)=c("EntrezGeneID","GeneSymbol","ENSG",paste(as.character(tumorInfo$patient),".FC",sep="") )
    return(fc.df)
}



## prepare count data and design matrix
myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
myContrast<-read.table(contrastFile,header=T) 
myCountTable<-readCountTable(myDesign,countDir)
head(myContrast)
## DESeq2 processing
if(args$pair)
{
	dds <- DESeqDataSetFromMatrix(countData = myCountTable[,c(-1,-2)], colData = myDesign, design = ~ patient + sampleType)
}else
{
	if(args$library)
	{
		dds <- DESeqDataSetFromMatrix(countData = myCountTable[,c(-1,-2)], colData = myDesign, design = ~ libType + sampleType)
	}else
	{
		dds <- DESeqDataSetFromMatrix(countData = myCountTable[,c(-1,-2)], colData = myDesign, design = ~ sampleType)
	}
}
## detect differential expressed genes
#register(MulticoreParam(8))
cat(sprintf("%s\tSetting up multicore (%d) ...\n", Sys.time(),args$nCPU))
register(MulticoreParam(args$nCPU))
dds <- DESeq(dds)
save.image(file=paste(outDir,"/DESeq2.RData",sep=""))

## QA and visulization
cat(sprintf("%s\tQC and visualization\n", Sys.time()))
vsd.blind <- varianceStabilizingTransformation(dds, blind=T)
vstMat.blind <- assay(vsd.blind)
qcAndViz(vsd.blind,vstMat.blind,myDesign,outDir,extra=paste0(".",sample.id))
save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
	
cat(sprintf("%s\tvst format for downstream...\n", Sys.time()))
vsd <- varianceStabilizingTransformation(dds, blind=F)
vstMat <- assay(vsd)
save.image(file=paste(outDir,"/DESeq2.RData",sep=""))

t.outdir <- outDir
runOneContrast <- function(l.ref,l.alt)
{
	#l.ref <- "norml"
	#l.alt <- "tumor"
	cat(sprintf("%s\tcontrast %s .vs. %s\n", Sys.time(),l.ref,l.alt))
	outDir <- paste0(t.outdir,"/",l.ref,"-",l.alt)
	print(paste0("making dir ",outDir))
	dir.create(outDir,recursive = T,showWarnings = F)

	res <- results(dds,contrast=c("sampleType",l.alt,l.ref))
	res
	#res <- results(dds)
	resOrdered <- res[order(res$padj),]
	resSig <- subset(resOrdered, padj < 0.1)
	resSigStrict <- subset(resSig, padj < opt.qval & abs(log2FoldChange) > 1)
	summary(resSigStrict)
	## output; for downstream analysis, blind=F
	#cat(sprintf("%s\tDE genes visualization\n", Sys.time()))
	## cluster
	#plotCluster(resSigStrict,vsd,myDesign,outDir)
	save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
	## per-sample fc
	fc.df<-data.frame()
	if(args$pair)
	{
	    fc.df<-makePersampleFC(vstMat,myDesign)
	    write.table(fc.df, file=paste(outDir,"/DESeq2.FCPerSample.txt",sep=""),sep="\t",quote=F,row.names=F)
	}
	cat(sprintf("%s\toutput report\n", Sys.time()))
	## DESeq2 report
	#des2Report <- outputDEGene(resSig,resSigStrict,dds,vstMat,fc.df,outDir)
	des2Report <- outputDEGene(resOrdered,resSigStrict,dds,vstMat,fc.df,outDir)
	cat(sprintf("%s\t---- des2Report done\n", Sys.time()))
	save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
	# GAGE report
	exp.fc <- res$log2FoldChange
	names(exp.fc) <- rownames(res)
	gageReport <- gageAnalysis(exp.fc,outDir)
	cat(sprintf("%s\t---- gageReport done\n", Sys.time()))
	save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
	##  GO and KEGG enrichment by hyperG test
	selectedIDs <- rownames(resSigStrict)
	universeIDs <- rownames(dds)
	goReport <- hyperGEA("GOHyperGParams",selectedIDs,universeIDs,"GO_BP",ontology = "BP",conditional = TRUE,outDir)
	cat(sprintf("%s\t---- goReport done\n", Sys.time()))
	KEGGReport <- hyperGEA("KEGGHyperGParams",selectedIDs,universeIDs,"KEGG",outDir)
	cat(sprintf("%s\t---- KEGGReport done\n", Sys.time()))
	save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
	## assembly all report
	indexPage <- HTMLReport(shortName = "indexRNASeq", title = "Analysis of RnaSeqData", reportDirectory="reports", basePath=outDir)
	publish(hwrite("Differential Expressed genes detected:", heading=4), indexPage)
	publish(Link(list(des2Report), report = indexPage), indexPage)
	publish(hwrite("Gene set enrichment analysis:", heading=4), indexPage)
	publish(Link(list(goReport, KEGGReport, gageReport), report = indexPage), indexPage)
	finish(indexPage)

	save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
}
### TEST
#runOneContrast("normal","tumor")
apply(myContrast,1,function(x){ runOneContrast(x[1],x[2]) } )

save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
cat(sprintf("%s\trun successfully!\n", Sys.time()))
#
