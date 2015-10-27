#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-c", "--countDir", type="character", required=TRUE, help="HTSeqGenie output directory")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-p", "--pair", action="store_true", default=FALSE, help="if specified samples are paired ")
#parser$add_argument("-n", "--add_numbers", action="store_true", default=FALSE, help="Print line number at the beginning of each line [default]")
#parser$add_argument("infile", nargs=1, help="Infile")
args <- parser$parse_args()
designFile <- args$designFile
countDir <- args$countDir
outDir <- args$outDir

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
readGMT<-function(file) 
{
    f <- readLines(file)
    lst = sapply(f, function(x) unlist(strsplit(x, "\t", fixed = TRUE)))
    names(lst) = sapply(lst, function(x) x[1])
    gSet = lapply(lst, function(x) x[-(1:2)])
    gLink = unlist(lapply(lst, function(x) x[2]))
    list(gSet=gSet,gLink=gLink)
}

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

# qcAndViz(dds,res,resSigStrict,myDesign,outDir)
qcAndViz<-function(x,res,resSigStrict,designM,outDir)
{
    ## rlog and VST; for QC, bind=T
    rld <- rlog(x, blind=T)
    vsd <- varianceStabilizingTransformation(x, blind=T)
    rlogMat <- assay(rld)
    vstMat <- assay(vsd)

    ## visualization
    color.map <- function(sampleType) { if (sampleType=="normal" || sampleType=="control") "#0000FF" else "#FF0000" }
    patientcolors <- unlist(lapply(designM$sampleType,color.map))
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

    pdf(paste(outDir,"/DESeq2.plot.QA.pdf",sep=""))
    ## MA plot
    plotMA(res, main="DESeq2", ylim=c(-2,2))
    ## dispersion plot
    plotDispEsts(x)
    ## gene-sample heatmap
    select <- order(rowMeans(counts(x,normalized=TRUE)),decreasing=TRUE)[1:50]
    dat.plot <-vstMat[select,]
    rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
    heatmap.2(dat.plot, Rowv = T, Colv = T, scale="row", dendrogram="both", trace="none", margin=c(10, 6), main="Most Highly Expressed genes",cexRow=0.5)    
    ## PCA
    if(args$pair)
    { 
    	plotPCA(vsd, intgroup=c("sampleType", "patient")) 
    }else
    {
    	plotPCA(vsd, intgroup=c("sampleType")) 
    }
    dev.off()

    ## heatmap by most variable genes
    rowVar <- apply(vstMat,1,var)
    select <- order(rowVar,decreasing = T)[1:50]
    dat.plot <-vstMat[select,] 
    rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
    png(paste(outDir,"/DESeq2.QC.Cluster.Var.png",sep=""),width=1200,height=1200)
    heatmap.2(dat.plot, Rowv = T, Colv = T, scale="row", dendrogram="both", trace="none", margin=c(10, 6), main="Most variable genes",cexRow=0.8)
    #heatmap.2(dat.plot, col = hmcol, Rowv = FALSE, Colv = T, scale="none", dendrogram="column", trace="none", margin=c(10, 6), main="Most variable genes")
    dev.off()

    ## sample distance
    distsRL <- dist(t(vstMat))
    mat <- as.matrix(distsRL)
    rownames(mat) <- colnames(mat) <- with(colData(x), paste(sampleType, patient, sep=" : "))
    png(paste(outDir,"/DESeq2.QC.SampleDist.png",sep=""),width=1200,height=1200)
    heatmap.2(mat, trace="none", margin=c(13, 13), main="Sample Distance")
    dev.off()
    #heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13), main="Sample Distance")
    #data <- plotPCA(vsd, intgroup=c("sampleType", "patient"), returnData=TRUE)
    #percentVar <- round(100 * attr(data, "percentVar"))
    #ggplot(data, aes(PC1, PC2, color=sampleType, shape=patient)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2:",percentVar[2],"% variance"))
    #dev.off()

    pdf(paste(outDir,"/DESeq2.plot.count.pdf",sep=""))
    ##
    select <- rownames(resSigStrict)[1:30]
    select.f <- !is.na(select)
    select <- select[select.f]
    print(select)
    if(length(select)>0)
    {
	for(i in 1:length(select))
	{
	    #plotCounts(dds, gene=select[i], intgroup="sampleType", returnData=F, main=paste("Gene: ",unlist(as.list(org.Hs.egSYMBOL)[select[i]]),"(",select[i],"), padj=",resSigStrict[select[i],"padj"],sep=""))
	    gg <- plotCounts(x, gene=select[i], intgroup="sampleType", transform=T, returnData=T, main=paste("Gene: ",entrezToXXX(select[i]),"(",select[i],"), padj=",resSigStrict[select[i],"padj"],sep=""))
	    gg$patient <- designM[rownames(gg),"patient"]
	    gg2<-cast(gg,patient~sampleType,value="count")
	    #beeswarm(Exp ~ Fusion, data=inTable, add = T,pch=20)
	    stripchart(count ~ sampleType, data = gg, vertical = TRUE, method = "overplot", jitter=0.01,  pch = 21, col = "maroon", bg = "bisque", main=paste("Gene: ",entrezToXXX(select[i]),"(",select[i],"), padj=",resSigStrict[select[i],"padj"],sep=""), ylab="log2(normalized count)")
	    segments(1,gg2[,2],2,gg2[,3],col="maroon")
	    #ggplot(d, aes(x=sampleType, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(25,100,400))
	}
    }
    dev.off()
    #cat("qcAndViz: Good\n")
}

plotCluster<-function(resSigStrict,vsd,designM,outDir)
{
    color.map <- function(sampleType) { if (sampleType=="normal" || sampleType=="control") "#0000FF" else "#FF0000" }
    patientcolors <- unlist(lapply(designM$sampleType,color.map))
    
    pdf(paste(outDir,"/DESeq2.plot.DECluster.pdf",sep=""))
    
    select <- rownames(resSigStrict)[1:30]
    select.f <- !is.na(select)
    select <- select[select.f]
    if(length(select)>0)
    {
    	heatmap.2(assay(vsd)[select,], col = greenred(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", trace="none", margin=c(10, 6), main="Most Differentially Expressed genes")
    }
    
    select <- rownames(resSigStrict)
    select.f <- !is.na(select)
    select <- select[select.f]
    if(length(select)>0)
    {
    	heatmap.2(assay(vsd)[select,], col = greenred(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", trace="none", margin=c(10, 6), main="All Differentially Expressed genes",labRow=F)
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

gageAnalysis<-function(exp.fc,outDir)
{
    data(kegg.gs)
    c2.cp.gs=readGMT("/Share/BP/zhenglt/00.database/MSigDB/msigdb_v4.0_files_to_download_locally/msigdb_v4.0_GMTs/c2.cp.v4.0.entrez.gmt")
    #c2.cp.gs=readList("/Share/BP/zhenglt/00.database/msigdb/msigdb_v4.0_files_to_download_locally/msigdb_v4.0_GMTs/c2.cp.v4.0.entrez.gmt")
    #### gage using kegg.gs
    gagePathway<-function(exp.fc,outDir,gs,gs.name,gl=NULL)
    {
	# 1 direction
	fc.p <- gage(exp.fc, gsets = gs, ref = NULL, samp = NULL)
	fc.p$greater=as.data.frame(fc.p$greater)
	fc.p$less=as.data.frame(fc.p$less)
	sel.h <- fc.p$greater[, "q.val"] < 0.25 & !is.na(fc.p$greater[, "q.val"])
	path.ids.h <- rownames(fc.p$greater)[sel.h]
	sel.l <- fc.p$less[, "q.val"] < 0.25 & !is.na(fc.p$less[,"q.val"])
	path.ids.l <- rownames(fc.p$less)[sel.l]
	# 2 direction
	fc.p.2d <- gage(exp.fc, gsets = gs, ref = NULL, samp = NULL, same.dir = F)
	fc.p.2d$greater=as.data.frame(fc.p.2d$greater)
	fc.p.2d$less=as.data.frame(fc.p.2d$less)
	sel.2d <- fc.p.2d$greater[, "q.val"] < 0.25 & !is.na(fc.p.2d$greater[, "q.val"])
	path.ids.2d <- rownames(fc.p.2d$greater)[sel.2d]
	# all path ids
	path.ids <- c(path.ids.h, path.ids.l, path.ids.2d)
	if(sum(regexpr("hsa[0-9]{5}",names(gs),perl=T)==rep(1,length(names(gs))))==length(names(gs)))
	{
		path.ids <- substr(path.ids, 1, 8)
	}
	## output table
    	out1<-fc.p$greater[sel.h,]
    	out1<-data.frame(pathway=rownames(out1),out1)
    	out2<-fc.p$less[sel.l,]
    	out2<-data.frame(pathway=rownames(out2),out2)
	out1d<-rbind(out1,out2)
    	out2d<-fc.p.2d$greater[sel.2d,]
    	out2d<-data.frame(pathway=rownames(out2d),out2d)
    	write.table(out1d, file=paste(outDir,"/DESeq2.gage.",gs.name,".1d.txt",sep=""),sep="\t",quote=F,row.names=F)
    	write.table(out2d, file=paste(outDir,"/DESeq2.gage.",gs.name,".2d.txt",sep=""),sep="\t",quote=F,row.names=F)
	####
	if(gs.name=="KEGG")
	{
       	    # pathview, only for "KEGG"
	    picDir <- paste(outDir,"/reports/figure_gage",sep="")
	    dir.create(picDir, showWarnings = FALSE)
	    out.suffix="gage.DESeq2"
	    oriDir <- getwd()
	    setwd(picDir)
	    ## OUT/reports/figure_gage/hsa04740.gage.DESeq2.png
	    pv.out.list <- sapply(path.ids, function(pid) pathview(kegg.dir="/Share/BP/zhenglt/00.database/kegg/pathview",gene.data =  exp.fc, pathway.id = pid, species = "hsa", out.suffix=out.suffix))
	    setwd(oriDir)
	    if(nrow(out1d) > 0)
	    {
		imagename <- paste0("figure_gage/",substr(out1d$pathway, 1, 8),".",out.suffix,".png")
		out1d$Image <- hwriteImage(imagename, link = imagename, table = FALSE, height=50, width=50)
	    }
	    if(nrow(out2d) > 0)
	    {
		imagename <- paste0("figure_gage/",substr(out2d$pathway, 1, 8),".",out.suffix,".png")
		out2d$Image <- hwriteImage(imagename, link = imagename, table = FALSE, height=50, width=50)
	    }
	}
	#### Report writting
	addGSetLink <- function(object, ...)
    	{
		if(!is.null(gl))
		{
    			object$pathway <- hwrite(as.character(object$pathway), link=gl[as.character(object$pathway)], table = FALSE)
		}
		return(object)
    	}
    	gageReport <- HTMLReport(shortName = paste("gage_analysis_rnaseq_",gs.name,"",sep=""), title = paste("GAGE analysis of RnaSeqData (", gs.name, ")",sep=""), reportDirectory = "reports", basePath=outDir)
    	publish(hwrite("One direction (up or down)", heading=2), gageReport)
	if(nrow(out1d) > 0)
    	{
		publish(out1d, gageReport, reportDir="./reports", .modifyDF=list(addGSetLink))
    	}
    	publish(hwrite("Two direction (up and down)", heading=2), gageReport)
    	if(nrow(out2d) > 0)
    	{
		publish(out2d, gageReport, reportDir="./reports", .modifyDF=list(addGSetLink))
    	}
    	finish(gageReport)
	list(report=gageReport)
    }
    gage.kegg=gagePathway(exp.fc,outDir,kegg.gs,"KEGG")
    gage.Canonical=gagePathway(exp.fc,outDir,c2.cp.gs$gSet,"Canonical",c2.cp.gs$gLink)
    ## return report obj for "index" web
    list(gage.kegg$report, gage.Canonical$report)
}

outputDEGene<-function(res,resSigStrict,dds,vstMat,fc.df,outDir)
{
    ## txt report
    resSig <- subset(res, padj < 0.1)
    resSigStrict <- subset(resSig, padj < 0.01 & abs(log2FoldChange) > 1)

    out_res <- as.data.frame(res)
    out_res$EntrezGeneID <- rownames(out_res)
    out_res$GeneSymbol <- entrezToXXX(out_res$EntrezGeneID)
    out_res$ENSG <- entrezToXXX(out_res$EntrezGeneID,"ENSG")
    
    if(nrow(fc.df)>0)
    {
    	out_res <- cbind(out_res, as.data.frame(vstMat)[rownames(out_res),], fc.df[rownames(out_res),c(-1,-2,-3)])
    }else
    {
    	out_res <- cbind(out_res, as.data.frame(vstMat)[rownames(out_res),])
    }
    
    out_resSig <- subset(out_res, padj < 0.1 )
    out_resSigStrict <- subset(out_resSig, padj < 0.01 & abs(log2FoldChange) > 1 )
    write.table(out_res, file=paste(outDir,"/DESeq2.txt",sep=""),sep="\t",quote=F,row.names=F)
    write.table(out_resSig, file=paste(outDir,"/DESeq2.sig.txt",sep=""),sep="\t",quote=F,row.names=F)
    write.table(out_resSigStrict, file=paste(outDir,"/DESeq2.sig.strict.txt",sep=""),sep="\t",quote=F,row.names=F)
    ## HTML report
    des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2', title = 'RNA-seq analysis of differential expression using DESeq2', reportDirectory="reports", basePath=outDir)
    if(nrow(resSigStrict)>0)
    {
	    publish(resSigStrict,des2Report, pvalueCutoff=0.01,lfc=1, n=Inf, annotation.db="org.Hs.eg.db", DataSet=dds, factor = colData(dds)$sampleType, reportDir="./reports")
    }
    finish(des2Report)
    des2Report
}

hyperGEA<-function(x,selectedIDs,universeIDs,reportName, ...)
{
    addGOIDLink <- function(object, ...)
    {
    	object$GOID <- hwrite(as.character(object$GOID), link = paste0("http://amigo.geneontology.org/amigo/term/", as.character(object$GOID)), table = FALSE)
	return(object)
    }
    addKEGGIDLink <- function(object, ...)
    {
    	#http://www.genome.jp/kegg/pathway/hsa/hsa05219.html
    	object$KEGGID <- hwrite(as.character(object$KEGGID), link = paste0("http://www.genome.jp/kegg/pathway/hsa/hsa", as.character(object$KEGGID),".html"), table = FALSE)
	return(object)
    }

    aReport <- HTMLReport(shortName = paste("Hypergeometric Tests for Gene Set (",reportName,")",sep=""), title = paste("Hypergeometric Tests for Gene Set (",reportName,")",sep=""), reportDirectory="reports", basePath=outDir)

    if(length(selectedIDs)>0)
    {
	aParams <- new(x, geneIds = selectedIDs, universeGeneIds = universeIDs, annotation ="org.Hs.eg", pvalueCutoff = 0.01, testDirection = "over", ... )
	aResults <- hyperGTest(aParams)
	aResultsSummary <- summary(aResults)

	if(nrow(aResultsSummary)>0)
	{
		if(x=="GOHyperGParams")
		{
		    #http://amigo.geneontology.org/amigo/term/GO:0000070
		    names(aResultsSummary)[1]<-"GOID"
		    publish(aResultsSummary, aReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg", pvalueCutoff= 0.01, .modifyDF=list(addGOIDLink))
		}else if(x=="KEGGHyperGParams")
		{
		    publish(aResultsSummary, aReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg", pvalueCutoff= 0.01, .modifyDF=list(addKEGGIDLink))
		}
	}
    }
    finish(aReport)
    aReport
   
    # not work still
    #pfamParams <- new("PFAMHyperGParams", geneIds= selectedIDs, universeGeneIds=universeIDs, annotation="org.Hs.eg", pvalueCutoff= 0.01, testDirection="over")
    #PFAMResults <- hyperGTest(pfamParams)
    #PFAMReport <- HTMLReport(shortName = 'pfam_analysis_rnaseq', title = "PFAM analysis of RnaSeqData", reportDirectory = paste(outDir,"/reports",sep=""))
    #publish(PFAMResults, PFAMReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg",categorySize=5)
    #finish(PFAMReport)

}

## prepare count data and design matrix
myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
myCountTable<-readCountTable(myDesign,countDir)
## DESeq2 processing
if(args$pair)
{
	dds <- DESeqDataSetFromMatrix(countData = myCountTable[,c(-1,-2)], colData = myDesign, design = ~ patient + sampleType)
}else
{
	dds <- DESeqDataSetFromMatrix(countData = myCountTable[,c(-1,-2)], colData = myDesign, design = ~ sampleType)
}
## detect differential expressed genes
register(MulticoreParam(4))
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.1)
resSigStrict <- subset(resSig, padj < 0.01 & abs(log2FoldChange) > 1)
summary(resSigStrict)
## QA and visulization
qcAndViz(dds,res,resSigStrict,myDesign,outDir)
## output; for downstream analysis, blind=F
rld <- rlog(dds, blind=F)
vsd <- varianceStabilizingTransformation(dds, blind=F)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
## cluster
plotCluster(resSigStrict,vsd,myDesign,outDir)
save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
## per-sample fc
fc.df<-data.frame()
if(args$pair)
{
    fc.df<-makePersampleFC(vstMat,myDesign)
    write.table(fc.df, file=paste(outDir,"/DESeq2.FCPerSample.txt",sep=""),sep="\t",quote=F,row.names=F)
}
## DESeq2 report
#des2Report <- outputDEGene(resSig,resSigStrict,dds,vstMat,fc.df,outDir)
des2Report <- outputDEGene(resOrdered,resSigStrict,dds,vstMat,fc.df,outDir)
cat("================ des2Report done ================\n")
save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
# GAGE report
exp.fc <- res$log2FoldChange
names(exp.fc) <- rownames(res)
gageReport <- gageAnalysis(exp.fc,outDir)
cat("================ gageReport done ================\n")
save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
##  GO and KEGG enrichment by hyperG test
selectedIDs <- rownames(resSigStrict)
universeIDs <- rownames(dds)
goReport <- hyperGEA("GOHyperGParams",selectedIDs,universeIDs,"GO_BP",ontology = "BP",conditional = TRUE)
cat("================ goReport done ================\n")
KEGGReport <- hyperGEA("KEGGHyperGParams",selectedIDs,universeIDs,"KEGG")
cat("================ KEGGReport done ================\n")
save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
## assembly all report
indexPage <- HTMLReport(shortName = "indexRNASeq", title = "Analysis of RnaSeqData", reportDirectory="reports", basePath=outDir)
publish(hwrite("Differential Expressed genes detected:", heading=4), indexPage)
publish(Link(list(des2Report), report = indexPage), indexPage)
publish(hwrite("Gene set enrichment analysis:", heading=4), indexPage)
publish(Link(list(goReport, KEGGReport, gageReport), report = indexPage), indexPage)
finish(indexPage)

save.image(file=paste(outDir,"/DESeq2.RData",sep=""))
