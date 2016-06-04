#!/usr/bin/env Rscript

#rdata.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/DESeq2.OUT.phase01-21.multiCore.postTCellMarkerFilter/P0205/DESeq2.RData"
#suppressPackageStartupMessages(library("R.utils"))
#lenv <- loadToEnv(rdata.file)
### make test data
#iTest <- c(grep("^PTC",colnames(lenv$myCountTable))[1:20],grep("^TTC",colnames(lenv$myCountTable))[1:20])
#myCountTable <- lenv$myCountTable[,c(1,2,iTest)]
#myCountTable[1:4,1:8]
#myDesign <- lenv$myDesign[colnames(lenv$myCountTable)[iTest],]
#save(myCountTable,myDesign,file="test.P0205.RData")

args <- commandArgs(T)
if(length(args)<5)
{
    cat(sprintf("run.scde.R <n.cores> <output.prefix> <rdata.file> <contrast> <cellType.color>\n"))
    q()
}
arg.n.cores <- as.integer(args[1])
output.prefix <- args[2]
rdata.file <- args[3]
arg.contrast <- args[4]
cellType.color.file  <- args[5]

#arg.n.cores <- 12
#output.prefix <- "./test"
#rdata.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/test.P0205.RData"
#myContrast <- c("PTC","TTC")

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
suppressPackageStartupMessages(library(scde))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require(boot))
setwd(dirname(output.prefix))

myContrast <- unlist(strsplit(arg.contrast,split="[,.]"))
suppressPackageStartupMessages(library("R.utils"))
lenv <- loadToEnv(rdata.file)
if(arg.contrast != "all") { 
    myDesign <- subset(lenv$myDesign,sampleType %in% myContrast) 
}else { 
    myDesign <- lenv$myDesign 
}
myCountTable <- cbind(lenv$myCountTable[,c(1,2)],lenv$myCountTable[,rownames(myDesign)])
output.rdata <- sprintf("%s.ediff.RData",output.prefix)
output.pdf <- sprintf("%s.pdf",output.prefix)
sampleTypeColor <- read.SampleTypeColor(cellType.color.file)

# clean up the dataset
cd <- myCountTable[,c(-1,-2)]
# omit genes that are never detected
cd <- cd[rowSums(cd)>0, ]
# omit cells with very poor coverage
cd <- cd[, colSums(cd)>1e4]

sg <- as.character(myDesign[colnames(cd),"sampleType"])
names(sg) <- colnames(cd)
sg <- as.factor(sg)

pdf(output.pdf,width = 10,height = 10)
par(mar=c(5,5,4,2),cex.lab=1.5)
# calculate models
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = arg.n.cores, threshold.segmentation = TRUE, save.crossfit.plots = TRUE, save.model.plots = TRUE, verbose = 1)
# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]
# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = TRUE)
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  sg, n.randomizations  =  100, n.cores  =  arg.n.cores, verbose  =  1)
# differentially expressed genes
ediff.out <- ediff[order(abs(ediff$cZ), decreasing = TRUE), ]
ediff.out.df <- data.frame(geneID=rownames(ediff.out),geneSymbol=entrezToXXX(rownames(ediff.out)))
ediff.out.df <- cbind(ediff.out.df,ediff.out)
ediff.out.df.sig <- subset(ediff.out.df,abs(cZ)>=1.96 & abs(mle)>1)
#ediff.out.df.sig <- subset(ediff.out.df,abs(cZ)>=1.5 & abs(mle)>1)
# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff.out.df, file = sprintf("%s.scde.ediff.out",output.prefix), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(ediff.out.df.sig, file = sprintf("%s.scde.ediff.sig.out",output.prefix), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

####### DESeq2, for comparison ###
#####suppressPackageStartupMessages(library("DESeq2"))
#####suppressPackageStartupMessages(library("BiocParallel"))
#####dds <- DESeqDataSetFromMatrix(countData = cd, colData = myDesign[colnames(cd),], design = ~ sampleType)
#####cat(sprintf("%s\tSetting up multicore (%d) ...\n", Sys.time(),arg.n.cores))
#####register(MulticoreParam(arg.n.cores))
#####dds <- DESeq(dds,parallel=TRUE)
#####res <- results(dds,contrast=c("sampleType",myContrast[2],myContrast[1]),parallel=TRUE)
#####resOrdered <- res[order(res$padj),]
#####DESeq2.out <- data.frame(geneID=rownames(resOrdered),geneSymbol=entrezToXXX(rownames(resOrdered)))
#####DESeq2.out <- cbind(DESeq2.out,resOrdered)
#####DESeq2.out.df <- DESeq2.out
#####DESeq2.out.df.sig <- subset(DESeq2.out,padj < 0.05 & abs(log2FoldChange)>1)
#####write.table(DESeq2.out.df, file = sprintf("%s.DESeq2.ediff.out",output.prefix), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#####write.table(DESeq2.out.df.sig, file = sprintf("%s.DESeq2.ediff.sig.out",output.prefix), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# get expression magntiude estimates
o.fpm <- scde.expression.magnitude(o.ifm, counts = cd)
# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm, magnitudes = log((10^o.prior$x)-1))
par(mfrow = c(1,1), mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1)
plot(c(), c(), xlim=range(o.prior$x), ylim=c(0,1), xlab="expression magnitude (log10)", ylab="drop-out probability")
invisible(apply(o.fail.curves[, grep(myContrast[1],colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y,col = "orange")))
invisible(apply(o.fail.curves[, grep(myContrast[2], colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y, col = "dodgerblue")))

dev.off()

# get self-fail probabilities (at a given observed count)
p.self.fail <- scde.failure.probability(models = o.ifm, counts = cd)
# simulate drop-outs
# note: using 10 sampling rounds for illustration here. ~500 or more should be used.
n.simulations <- 10; k <- 0.9;
cell.names <- colnames(cd); names(cell.names) <- cell.names;
dl <- mclapply(1:n.simulations,function(i) {
  scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
    x <- cd[,nam];
    # replace predicted drop outs with NAs
    x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
    x;
  }))
  rownames(scd1) <- rownames(cd); 
  # calculate correlation on the complete observation pairs
  cor(log10(scd1+1),use="pairwise.complete.obs");
}, mc.cores = arg.n.cores)
# calculate average distance across sampling rounds
direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))

# load boot package for the weighted correlation implementation
k <- 0.95;
reciprocal.dist <- as.dist(1 - do.call(rbind, mclapply(cell.names, function(nam1) {
  unlist(lapply(cell.names, function(nam2) {
    # reciprocal probabilities
    f1 <- scde.failure.probability(models = o.ifm[nam1,,drop = FALSE], magnitudes = o.fpm[, nam2])
    f2 <- scde.failure.probability(models = o.ifm[nam2,,drop = FALSE], magnitudes = o.fpm[, nam1])
    # weight factor
    pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
    boot::corr(log10(cbind(cd[, nam1], cd[, nam2])+1), w = pnf)
  }))
},mc.cores = arg.n.cores)), upper = FALSE)

# reclculate posteriors with the individual posterior modes 
jp <- scde.posteriors(models = o.ifm, cd, o.prior, return.individual.posterior.modes = TRUE, n.cores = arg.n.cores)
# find joint posterior modes for each gene - a measure of MLE of group-average expression
jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
p.mode.fail <- scde.failure.probability(models = o.ifm, magnitudes = jp$jp.modes)
# weight matrix
matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
# magnitude matrix (using individual posterior modes here)
mat <- log10(exp(jp$modes)+1);
# weighted distance
mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
  unlist(lapply(cell.names,function(nam2) {
    corr(cbind(mat[, nam1], mat[, nam2]), w = sqrt(sqrt(matw[, nam1]*matw[, nam2])))
  }))
}, mc.cores = arg.n.cores)), upper = FALSE)

pdf(sprintf("%s.dist.pdf",output.prefix),width = 10,height = 10)
par(mar=c(5,5,4,2),cex.main=1.8)
plot(hclust(direct.dist),cex.lab=1.8,cex=0.4*120/nrow(myDesign))
plot(hclust(reciprocal.dist),cex.lab=1.8,cex=0.4*120/nrow(myDesign))
plot(hclust(mode.fail.dist),cex.lab=1.8,cex=0.4*120/nrow(myDesign))
dev.off()

save.image(file = output.rdata)

