#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-c", "--countDir", type="character", required=TRUE, help="HTSeqGenie output directory")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-l", "--cyclebase", type="character", default="/Share/BP/zhenglt/02.pipeline/cancer/rna/data/cycleBase.human.RData", help="sample id [default %(default)s]")
parser$add_argument("-p", "--process", action="store_true", default=FALSE, help="if specified run scLVM analysis; otherwise  preprocess data only, the scLVM run shoud be done using python script [default %(default)s]")
args <- parser$parse_args()
designFile <- args$designFile
countDir <- args$countDir
out.dir <- args$outDir
sample.id <- args$sample
cyclebase.rdata <- args$cyclebase
do.scLVM <- args$process

#designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/sample.design.Tang.P1202.S126"
#countDir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/OUT/Tang.P1202"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/scLVM"
#sample.id <- "Tang"

dir.create(out.dir,recursive = T,showWarnings = F)

suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(statmod))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(HTSeqGenie))
suppressPackageStartupMessages(library(scLVM))
suppressPackageStartupMessages(library(rhdf5))

## function definition
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
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

myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
myCountTable<-readCountTable(myDesign,countDir)
#myContrast<-read.table(contrastFile,header=T) 

countData = myCountTable[,c(-1,-2)]
countData.sf <- estimateSizeFactorsForMatrix(countData)
#normalise read counts
countData.sfNormalized <- t( t(countData) / countData.sf )

#get technical noise
pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noSpikein.log.pdf"),width = 16,height = 8)
techNoiseLogFit = fitTechnicalNoise(countData.sfNormalized, fit_type = 'log', use_ERCC = FALSE, plot=TRUE) 
dev.off()

pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noSpikein.logvar.pdf"),width = 16,height = 8)
techNoiseLogVarFit = fitTechnicalNoise(countData.sfNormalized, fit_type = 'logvar', use_ERCC = FALSE, plot=TRUE) 
dev.off()

#call variable genes
pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noSpikein.log.variableGenes.pdf"),width = 16,height = 8)
is_hetLog = getVariableGenes(countData.sfNormalized, techNoiseLogFit$fit, method = "fit", plot=TRUE)
dev.off()
table(is_hetLog)

pdf(paste0(out.dir,"/",sample.id,".fitTechNoise.noSpikein.logvar.variableGenes.pdf"),width = 16,height = 8)
is_hetLogVar = getVariableGenes(countData.sfNormalized, techNoiseLogVarFit$fit, method = "fit", plot=TRUE)
dev.off()
table(is_hetLogVar)

#get cell cycle genes from GO 
## scLVM's bug: the name of the vecotr holding genes should be "ens_ids_cc"
#gid.cellCycle <- getEntrezFromGO("GO:0007049")
ens_ids_cc <- getEntrezFromGO("GO:0007049")
lname <- load(cyclebase.rdata)
cellcyclegenes_filter <- na.omit(match(ens_ids_cc,rownames(countData.sfNormalized)))
cellcyclegenes_filterCB <- na.omit(match(dataCB[1:600,"Entrez"],rownames(countData.sfNormalized)))

#rename a few variables
Y = log10(countData.sfNormalized+1) #normalised trandformed read counts
genes_het_bool = as.vector(is_hetLog) #variable genes
geneID = rownames(countData.sfNormalized) #gene IDs
tech_noise = as.vector(techNoiseLogFit$techNoiseLog) #technical noise
cell_names <- colnames(countData.sfNormalized)
write.table(geneID[genes_het_bool],file=paste0(out.dir,"/",sample.id,".var.geneID.txt"),quote = F,row.names = F)
write.table(cell_names,file=paste0(out.dir,"/",sample.id,".var.cellNames.txt"),quote = F,row.names = F)
write.table(rownames(countData.sfNormalized)[cellcyclegenes_filter],file=paste0(out.dir,"/",sample.id,".var.cellcyclegenesGO.txt"),quote = F,row.names = F)
write.table(rownames(countData.sfNormalized)[cellcyclegenes_filterCB],file=paste0(out.dir,"/",sample.id,".var.cellcyclegenesCB.txt"),quote = F,row.names = F)
h5save(cellcyclegenes_filter,cellcyclegenes_filterCB,geneID,cell_names,genes_het_bool,tech_noise,Y,file=paste0(out.dir,"/",sample.id,".LogCounts.h5f"))
###rm(list())
## hdf5 and R's data.frame have opposite dimenshion direction
Y = t(Y) 

if(do.scLVM)
{
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

