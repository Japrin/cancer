#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-c", "--countDir", type="character", required=TRUE, help="HTSeqGenie output directory")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-f", "--fdr", type="double", default="0.001", help="FDR threshold  [default %(default)s]")
parser$add_argument("-l", "--cyclebase", type="character", default="/Share/BP/zhenglt/02.pipeline/cancer/rna/data/cycleBase.human.RData", help="sample id [default %(default)s]")
parser$add_argument("-e", "--ercc", action="store_true", default=FALSE, help="if specified normalize endo-gene using ERCCs' size factor; otherwise normalize endo-gene using endo-genes' size factor  [default %(default)s]")
parser$add_argument("-n", "--ignoreERCC", action="store_true", default=FALSE, help="whether ignore ERCC data even it's available [default %(default)s]")
parser$add_argument("-x", "--excludeSamples", type="character", help="file contain sample list or comma sperated string, specify the samples to be excluded")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-p", "--process", action="store_true", default=FALSE, help="if specified run scLVM analysis; otherwise  preprocess data only, the scLVM run shoud be done using python script [default %(default)s]")
parser$add_argument("-t", "--tSNE", action="store_true", default=FALSE, help="if specified run tSNE analysis; otherwise  preprocess data only [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", type="character", help="cellTypeColorFile")
args <- parser$parse_args()
designFile <- args$designFile
countDir <- args$countDir
out.dir <- args$outDir
sample.id <- args$sample
fit.fdr <- args$fdr
cyclebase.rdata <- args$cyclebase
use.ERCC.sf <- args$ercc
do.scLVM <- args$process
mode.verbose <- args$verbose
ignore.ERCC <- args$ignoreERCC
exclude.samples <- args$excludeSamples
run.tSNE <- args$tSNE
cellTypeColorFile <- args$cellTypeColorFile

print(args)

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

getEntrezFromGO <- function(term)
{
	if (require(org.Hs.eg.db)) 
	{
		xxGO <- AnnotationDbi::as.list(org.Hs.egGO2EG)
	}
	else
	{
		stop("Install org.Hs.eg.db package for retrieving gene lists from GO")
	}
   	cell.cycleEG <- unlist(xxGO[term])
	cell.cycleEG
}
# plot size factor distribution
plotSizeFactorDist <- function(sf,out.prefix,sample.id.toHighlight=NULL)
{
    pdf(sprintf("%s.%s",out.prefix,"sizeFactor.pdf"),width = 8,height = 6)
    par(cex.lab=1.8,cex.axis=1.5,mar=c(5,6,4,2)+0.1)
    cdata <- sort(sf)
    ccol <- rep("darkblue",length(sf))
    names(ccol) <- names(cdata)
    if(!is.null(sample.id.toHighlight)) {
        ccol[sample.id.toHighlight] <- "red"
    }else{
        ccol[cdata < 0.4] <- "red"
    }
    barplot(cdata,col=ccol,border=NA,xlab="cell index",ylab="size factor",xaxt="n")
    box(which = "inner",lwd=4)
    par(new=TRUE, oma=c(1,8,1,1),mar=c(5,4,1,2)+0.1)
    layout(matrix(4:1,2))
    b.midpoint <- barplot(cdata[1:10],col=ccol[1:10],border=NA,xlab="",ylab="",xaxt="n")
    text(b.midpoint,y=-0.04, srt = 45, adj = 1, labels = names(cdata)[1:10], xpd = TRUE,cex=1.0)
    box(which = "figure")
    dev.off()
    pdf(sprintf("%s.%s",out.prefix,"sizeFactor.dist.pdf"),width = 8,height = 8)
    hist(sf,breaks = 30,freq = F,col = "gray",xlab="Size Factor",main=sprintf("%s",sample.id))
    lines(density(sf),col="orange",lwd=2)
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

plotSizeFactorDist(obj.scdn@size.factor.endo,sprintf("%s/%s.endo",out.dir,sample.id),sample.id.toHighlight=NULL)
out.size.factor.df <- data.frame(cellName=names(obj.scdn@size.factor.endo),szieFactor=obj.scdn@size.factor.endo)
write.table(out.size.factor.df,sprintf("%s/%s.endo.sizeFactor.txt",out.dir,sample.id),sep = "\t",row.names = F,col.names = T,quote = F)
print(quantile(obj.scdn@size.factor.endo,probs=c(0.001,0.002,0.005,0.01,0.05,0.1,0.5,0.9,0.95,0.99,0.995,0.998,0.999)))
if(obj.scdn@withERCC) {
    plotSizeFactorDist(obj.scdn@size.factor.ERCC,sprintf("%s/%s.ERCC",out.dir,sample.id))
    out.size.factor.df <- data.frame(cellName=names(obj.scdn@size.factor.ERCC),szieFactor=obj.scdn@size.factor.ERCC)
    write.table(out.size.factor.df,sprintf("%s/%s.ERCC.sizeFactor.txt",out.dir,sample.id),sep = "\t",row.names = F,col.names = T,quote = F)
    print(quantile(obj.scdn@size.factor.ERCC,probs=c(0.001,0.002,0.005,0.01,0.05,0.1,0.5,0.9,0.95,0.99,0.995,0.998,0.999)))
}

#get technical noise
pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.ERCC.counts.pdf"),width = 8,height = 8)
techNoiseERCCCounts = SCDenoise.fitTechnicalNoise(obj.scdn, fit_type = 'counts',plot=TRUE)
dev.off()
###q()
if(mode.verbose)
{
    pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noERCC.counts.pdf"),width = 8,height = 8)
    techNoiseEndoCounts = SCDenoise.fitTechnicalNoise(obj.scdn, fit_type = 'counts',use_ERCC=FALSE,plot=TRUE)
    dev.off()

    pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noERCC.log.pdf"),width = 8,height = 8)
    techNoiseLogFit = SCDenoise.fitTechnicalNoise(obj.scdn, fit_type = 'log', use_ERCC = FALSE, plot=TRUE) 
    dev.off()

    pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noERCC.logvar.pdf"),width = 8,height = 8)
    techNoiseLogVarFit = SCDenoise.fitTechnicalNoise(obj.scdn, fit_type = 'logvar', use_ERCC = FALSE, plot=TRUE) 
    dev.off()

    #call variable genes
    pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noERCC.counts.variableGenes.pdf"),width = 8,height = 8)
    vg.EndoCounts = SCDenoise.getVariableGenes(obj.scdn, techNoiseEndoCounts$fit, method = "fit", plot=TRUE)
    is_het_EndoCounts <- vg.EndoCounts[["is_het"]]
    dev.off()
    table(is_het_EndoCounts)
    vg.EndoCounts.df <- cbind(geneID=names(vg.EndoCounts[["is_het"]]),as.data.frame(vg.EndoCounts))
    vg.EndoCounts.df$geneSymbol <- entrezToXXX(vg.EndoCounts.df$geneID)
    vg.EndoCounts.df <- vg.EndoCounts.df[order(-vg.EndoCounts.df$residual),]
    write.table(vg.EndoCounts.df,sprintf("%s/%s.fitTechNoise.noERCC.counts.variableGenes.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)

    pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noERCC.log.variableGenes.pdf"),width = 8,height = 8)
    vg.Log = SCDenoise.getVariableGenes(obj.scdn, techNoiseLogFit$fit, method = "fit", plot=TRUE)
    is_hetLog <- vg.Log[["is_het"]]
    dev.off()
    table(is_hetLog)

    pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noERCC.logvar.variableGenes.pdf"),width = 8,height = 8)
    vg.LogVar = SCDenoise.getVariableGenes(obj.scdn, techNoiseLogVarFit$fit, method = "fit", plot=TRUE)
    is_hetLogVar <- vg.LogVar[["is_het"]]
    dev.off()
    table(is_hetLogVar)
   
    if(obj.scdn@withERCC)
    { 
        for(ff in unique(c(0.1,0.05,0.01,0.005,0.001,fit.fdr)))
        {
            pdf(sprintf("%s/%s.fitTechNoise.ERCC.counts.fdr%s.variableGenes.pdf",out.dir,sample.id,ff),width = 8,height = 8)
            vg.ERCCCounts = SCDenoise.getVariableGenes(obj.scdn, techNoiseERCCCounts$fit, method = "fdr", threshold = ff, fit_type="counts", plot=TRUE,fitB = techNoiseEndoCounts$fit)
            is_het <- vg.ERCCCounts[["is_het"]]
            dev.off()
            table(is_het)

            pdf(sprintf("%s/%s.fitTechNoise.venn.fdr%s.pdf",out.dir,sample.id,ff),width = 8,height = 8)
            venn.res <- venn(list(ERCC=names(is_het)[is_het],Endo=names(is_het_EndoCounts)[is_het_EndoCounts],log=names(is_hetLog)[is_hetLog],logVar=names(is_hetLogVar)[is_hetLogVar]))
            dev.off()
            cat(sprintf("%s ERCC .vs. Endo (fitting by counts), fdr: %s\n",sample.id,ff))
            table.ERCC.Endo <- table(ERCC=is_het,Endo=is_het_EndoCounts)
            print(addmargins(table.ERCC.Endo))
            print(round(sweep(addmargins(table.ERCC.Endo, 1, list(All = sum)), 2, apply(table.ERCC.Endo, 2, sum)/100, "/"), 2))
        }
    }else
    {
        pdf(sprintf("%s/%s.fitTechNoise.venn.noERCC.pdf",out.dir,sample.id),width = 8,height = 8)
        venn.res <- venn(list(Endo=names(is_het_EndoCounts)[is_het_EndoCounts],log=names(is_hetLog)[is_hetLog],logVar=names(is_hetLogVar)[is_hetLogVar]))
        dev.off()
    }
}

pdf(sprintf("%s/%s.fitTechNoise.ERCC.counts.variableGenes.pdf",out.dir,sample.id),width = 8,height = 8)
vg.ERCCCounts = SCDenoise.getVariableGenes(obj.scdn, techNoiseERCCCounts$fit, method = "fdr", threshold = fit.fdr, fit_type="counts", plot=TRUE)
is_het <- vg.ERCCCounts[["is_het"]]
dev.off()
table(is_het)
vg.ERCCCounts.df <- cbind(geneID=names(vg.ERCCCounts[["is_het"]]),as.data.frame(vg.ERCCCounts))
vg.ERCCCounts.df$geneSymbol <- entrezToXXX(vg.ERCCCounts.df$geneID)
if("padjA" %in% names(vg.ERCCCounts.df)) { 
    vg.ERCCCounts.df <- vg.ERCCCounts.df[order(vg.ERCCCounts.df$padjA,-vg.ERCCCounts.df$residual),]
}else
{
    vg.ERCCCounts.df <- vg.ERCCCounts.df[order(-vg.ERCCCounts.df$residual),]
}
write.table(vg.ERCCCounts.df,sprintf("%s/%s.fitTechNoise.ERCC.counts.variableGenes.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)


#rename a few variables
Y = log10(obj.scdn@normalized.endo+1) #normalised trandformed read counts
geneID = rownames(obj.scdn@normalized.endo) #gene IDs
cell_names <- colnames(obj.scdn@normalized.endo)
genes_het_bool = as.vector(is_het) #variable genes
tech_noise = as.vector(techNoiseERCCCounts$techNoiseLog) #technical noise

#goResults <- runTopGOAnalysis(geneID[genes_het_bool], geneID)
#print(goResults[["MF"]])
#print(goResults[["CC"]])
#print(goResults[["BP"]])

#get cell cycle genes from GO 
## scLVM's bug: the name of the vecotr holding genes should be "ens_ids_cc"
#gid.cellCycle <- getEntrezFromGO("GO:0007049")
ens_ids_cc <- getEntrezFromGO("GO:0007049")
lname <- load(cyclebase.rdata)
cellcyclegenes_filter <- na.omit(match(ens_ids_cc,geneID))
cellcyclegenes_filterCB <- na.omit(match(dataCB[1:600,"Entrez"],geneID))

write.table(geneID[genes_het_bool],file=paste0(out.dir,"/",sample.id,".var.geneID.txt"),quote = F,row.names = F)
write.table(cell_names,file=paste0(out.dir,"/",sample.id,".var.cellNames.txt"),quote = F,row.names = F)
write.table(geneID[cellcyclegenes_filter],file=paste0(out.dir,"/",sample.id,".var.cellcyclegenesGO.txt"),quote = F,row.names = F)
write.table(geneID[cellcyclegenes_filterCB],file=paste0(out.dir,"/",sample.id,".var.cellcyclegenesCB.txt"),quote = F,row.names = F)
h5save(cellcyclegenes_filter,cellcyclegenes_filterCB,geneID,cell_names,genes_het_bool,tech_noise,Y,file=paste0(out.dir,"/",sample.id,".LogCounts.h5f"))

out.df <- data.frame(geneID=geneID[genes_het_bool])
out.df$geneName <- entrezToXXX(out.df$geneID)
out.df <- cbind(out.df, obj.scdn@normalized.endo[genes_het_bool,])
write.table(out.df,file=sprintf("%s/%s.het.countGeneData.sfNormalized",out.dir,sample.id),sep="\t",quote = F,row.names = F)

save.image(file = sprintf("%s/%s.RData",out.dir,sample.id))

if(run.tSNE){
    sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
    qcAndVizMat(Y[is_het,],myDesign,out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=100000,extra=paste0(".",sample.id),sfilter=NULL, gfilter=NULL)
}
###rm(list())
if(do.scLVM)
{
    ## hdf5 and R's data.frame have opposite dimenshion direction
    Y = t(Y) 
    #construct and initialize new scLVM object
    sclvm = new("scLVM")
    sclvm = init(sclvm,Y=Y,tech_noise = tech_noise)

    CellCycleARD = fitFactor(sclvm,geneSet = ens_ids_cc, k=20,use_ard = TRUE)

    png(paste0(out.dir,"/",sample.id,".scLVM.cell-cycle.png"),width = 800,height = 800)
    plot(seq(1, length(CellCycleARD$X_ard)), CellCycleARD$X_ard, xlab = '# Factor', ylab = 'Variance explained')
    title('Variance explained by latent factors')
    dev.off()

    CellCycle = fitFactor(sclvm,geneSet = ens_ids_cc,k=1)

    #Get cell-cycle factor
    Kcc = CellCycle$K
    Xcc = CellCycle$X

    png(paste0(out.dir,"/",sample.id,".scLVM.similarityMatrix.pdf"),width = 800,height = 800)
    #Plot inferred similarity matrix
    image(Kcc,xaxt = "n", yaxt = "n", xlab = 'cells', ylab = 'cells')
    title('Similarity matrix based on cell cycle')
    dev.off()

    idx_het = which(is_hetLog)

    # fit the model for variable genes
    sclvm = varianceDecomposition(sclvm, K=Kcc, idx = idx_het)

    # get variance components
    results_var = getVarianceComponents(sclvm)
    var_filtered = results_var$var[results_var$conv,] # filter out genes for which vd has not converged
    head(var_filtered)

    # get corrected expression levels
    Ycorr = getCorrectedExpression(sclvm)
    dim(Ycorr)

    var_mean = apply(var_filtered,2,mean)
    colors = c('Green','Blue','Gray')
    pdf(paste0(out.dir,"/",sample.id,".scLVM.varianceComponets.pie.pdf"),width = 8,height = 8)
    pie(var_mean, , col = colors)
    dev.off()

    #idx_lmm = idx_het[1:50]
    idx_lmm = idx_het

    # fit lmm without correction
    res_nocorr = LMM(sclvm, K = NULL,idx = idx_lmm,verbose=TRUE)

    # fit lmm with correction

    res_corr = LMM(sclvm, K = Kcc, idx = idx_lmm,verbose = TRUE)


    png(paste0(out.dir,"/",sample.id,".scLVM.correlation.01.png"),width = 800,height = 800)
    heatmap.2(res_nocorr$beta, Rowv = NULL, Colv = NULL, dendrogram = "none",
                    labCol = as.character(idx_lmm), labRow = as.character(idx_lmm),srtCol = 0, key=T,density.info = "none",
                          trace="none", breaks=seq.int(from = -0.6, to = 1.0, length.out = 13), main = 'Without Correction')
    dev.off()

    png(paste0(out.dir,"/",sample.id,".scLVM.correlation.02.png"),width = 800,height = 800)
    heatmap.2(res_corr$beta, Rowv = NULL, Colv = NULL, dendrogram = "none",
                    labCol = as.character(idx_lmm), labRow = as.character(idx_lmm),srtCol = 0, key=T,density.info = "none",
                          trace="none", breaks=seq.int(from = -0.6, to = 1.0, length.out = 13), main = 'With Correction')
    dev.off()


    Yhet = Y[,idx_het]
    #geneSymbols = getSymbols(colnames(Yhet))
    geneSymbols <- entrezToXXX(colnames(Yhet))


    gene_plot = "GATA3"
    idx_gene = which(geneSymbols==gene_plot)

    #PCA on corrected data
    pcaCorr = prcomp(Ycorr,2)
    d <- qplot(pcaCorr$x[,1], pcaCorr$x[,2],colour=Ycorr[,idx_gene], xlab = 'PC1', ylab = 'PC2')
    d + ggtitle('PCA corrected gene expression') + scale_color_continuous(name =gene_plot)
    ggsave(file=paste0(out.dir,"/",sample.id,".PCA.corrected.png"))
    #PCA on uncorrected data
    pca = prcomp(Yhet,2)
    d <- qplot(pca$x[,1], pca$x[,2],colour=Yhet[,idx_gene], xlab = 'PC1', ylab = 'PC2')
    d + ggtitle('PCA uncorrected gene expression') + scale_color_continuous(name =gene_plot)
    ggsave(file=paste0(out.dir,"/",sample.id,".PCA.uncorrected.png"))

    save.image(file = paste0(out.dir,"/",sample.id,".scLVM.RData"))
    ####
    ######## Fitting multiple latent factors
    ###### First, let's generate a new scLVM object and fit the cell-cycle factor:
    #####get cell cycle genes from GO 
    ####ens_ids_cc <- getEnsembl('GO:0007049')
    ####
    #####construct and initialize new scLVM object
    ####sclvmMult = new("scLVM")
    ####sclvmMult = init(sclvmMult,Y=Y,tech_noise = tech_noise)
    ####
    ####CellCycle = fitFactor(sclvmMult,geneSet = ens_ids_cc, k=1)
    ####
    #####Get cell-cycle factor
    ####Kcc = CellCycle$K
    ####Xcc = CellCycle$X
    ####
    ####### After having fit the cell-cycle factor, we can now fit the Th2 factor by conditioning on the the dominant cell-cycle factor. As we might expect interactions between cell cycle and differentiation, we can also fit an interaction term.
    ####
    #####Load Th2 genes
    ####Th2_genes = read.table(system.file("extdata","Th2_markers.txt",package = "scLVM"), as.is=TRUE)$V1
    ####
    #####get Th2 marker genes 
    ####gene_symbols = getSymbols(rownames(dataMouse))
    ####idx_Th2 <- na.omit(match(Th2_genes, gene_symbols))
    ####
    ####th2 = fitFactor(sclvmMult, idx = idx_Th2, XKnown = Xcc, k = 1, interaction=TRUE)
    ####KTh2 = th2$K
    ####Kint = th2$Kint

}

