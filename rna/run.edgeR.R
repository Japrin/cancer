#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input file")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, 
                    help="design matrix file, should be sorted according sampleType")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="outDir")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id")
parser$add_argument("-y", "--sampleTypeColorFile", type="character",default="", help="sampleTypeColorFile")
#parser$add_argument("-b", "--contrastFile", type="character", required=TRUE, help="contrast file")
#parser$add_argument("-p", "--pair", action="store_true", default=FALSE, help="if specified samples are paired ")
#parser$add_argument("-l", "--library", action="store_true", default=FALSE, help="if specified consider libraryType as a confounding factor")
#parser$add_argument("-q", "--qval", default=0.01, type="double", help="qval for differential genes [default %(default)s]")
#parser$add_argument("-n", "--nCPU", type="integer", default=8,  help="Number of cpu to use [default %(default)s]")
#parser$add_argument("-n", "--add_numbers", action="store_true", default=FALSE, help="Print line number at the beginning of each line [default]")
args <- parser$parse_args()
inFile <- args$inFile
designFile <- args$designFile
outDir <- args$outDir
sample.id <- args$sample
sampleTypeColorFile <- args$sampleTypeColorFile
#contrastFile <- args$contrastFile
#opt.qval <- args$qval
#contrastFile <- "contrast.list"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
library("edgeR")

inFile <- "/WPS1/zhenglt/work/proj_xy/integrated/quantification/P1118.count.tab.gz"
designFile <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/clusteringResult/lung.leaf.addMajor.addCyclone.addClonotype.CD8.P1118.txt"
outDir <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/PNT/P1118"
sample.id <- "P1118"
sampleTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/CellType.color"
sampleType <- "sampleType"

dir.create(outDir,recursive = T,showWarnings = F)

myDesign <- read.delim(designFile,header = T,sep = "\t",check.names=F,stringsAsFactors=F)
rownames(myDesign) <- myDesign$sample
in.table <- read.table(inFile,header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F)
g.GNAME <- in.table[,1]
in.table <- as.matrix(in.table[,-1])
s.f <- intersect(rownames(myDesign),colnames(in.table))
myDesign <- myDesign[s.f,]
in.table <- in.table[,s.f]
mType <- length(unique(myDesign[,sampleType]))
if(!is.null(sampleTypeColorFile) && file.exists(sampleTypeColorFile)){
    sampleTypeColor <- read.SampleTypeColor(sampleTypeColorFile)
    sampleTypeColor <- sampleTypeColor[ names(sampleTypeColor) %in% unique(myDesign[,sampleType]) ]
}else{
    sampleTypeColor <- structure(colorRampPalette(brewer.pal(mType,"Paired"))(mType),
                              names=unique(myDesign[,sampleType]))
}

y <- DGEList(counts = in.table,
             genes = data.frame(geneID=rownames(in.table),geneSymbol=g.GNAME),
             samples = myDesign)
f <- apply(y$counts,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
y <- y[f,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)

pdf(sprintf("%s/%s.edgeR.MDS.pdf",outDir,sample.id),width = 8,height = 8)
plotMDS(y,labels = NULL,pch=16,col=sampleTypeColor[y$samples[,sampleType]])
dev.off()
y$samples$location <- relevel(factor(y$samples$location),ref = "P")
design <- model.matrix(~majorCluster+location,data = y$samples)

#register(MulticoreParam(args$nCPU))
y <- estimateDisp(y, design, robust=TRUE)

pdf(sprintf("%s/%s.edgeR.BCV.pdf",outDir,sample.id),width = 8,height = 8)
plotBCV(y)
dev.off()

fit <- glmFit(y, design)
for(sCoef in c("locationT","locationN"))
    { 
    lrt <- glmLRT(fit,coef = sCoef)
    topTags(lrt,n = 10)
    summary(de <- decideTestsDGE(lrt))
    out.df.t <- lrt$table[order(lrt$table$PValue),]
    out.df <- data.frame(geneID=rownames(out.df.t),
                         geneSymbol=y$genes[match(rownames(out.df.t),y$genes$geneID),"geneSymbol"]) 
    out.df <- cbind(out.df,out.df.t)
    out.df$FDR <- p.adjust(out.df$PValue,method = "BH")
    write.table(subset(out.df,abs(logFC)>1 & FDR<0.05),row.names = F,
                file = sprintf("%s/%s.edgeR.%s.sig.txt",outDir,sample.id,sCoef),
                sep = "\t",quote = F)
}

q()

