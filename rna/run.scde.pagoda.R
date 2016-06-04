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
    cat(sprintf("run.scde.R <n.cores> <output.prefix> <rdata.file> <contrast> <cellType.color> [use.gname]\n"))
    q()
}
arg.n.cores <- as.integer(args[1])
output.prefix <- args[2]
rdata.file <- args[3]
arg.contrast <- args[4]
cellType.color.file  <- args[5]

use.gname <- TRUE
if(length(args)>5) { use.gname <- FALSE }

#arg.n.cores <- 12
#output.prefix <- "./test"
#rdata.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/test.P0205.RData"
#arg.contrast <- "PTC,TTC"
#cellType.color.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/CellType.color"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
suppressPackageStartupMessages(library(scde))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require(boot))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GO.db))
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
output.rdata <- sprintf("%s.pagoda.RData",output.prefix)
output.pdf <- sprintf("%s.pdf",output.prefix)
sampleTypeColor <- read.SampleTypeColor(cellType.color.file)

# clean up the dataset
cd <- myCountTable[,c(-1,-2)]
# omit genes that are never detected
cd <- cd[rowSums(cd)>0, ]
# omit cells with very poor coverage
cd <- cd[, colSums(cd)>1e4]
if(use.gname)
{
    gname <- entrezToXXX(rownames(cd))
    cd <- cd[!is.na(gname),]
    rownames(cd) <- entrezToXXX(rownames(cd))
}
sg <- as.character(myDesign[colnames(cd),"sampleType"])
names(sg) <- colnames(cd)
sg <- as.factor(sg)

############### Pathway and Gene Set Overdispersion Analysis

pdf(sprintf("%s.pagoda.pdf",output.prefix),width = 20,height = 10)
par(mar=c(5,5,4,2),cex.lab=1.5)
# 
knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = arg.n.cores, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = arg.n.cores, plot = TRUE)
# list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

### make gene set from go
if(use.gname)
{
    ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
    rids <- names(ids)
    names(rids) <- ids
}else
{
    ids <- rownames(cd)
    rids <- ids
    names(rids) <- as.character(ids)
}
# list all the ids per GO category
go.env <- eapply(org.Hs.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
# omit categories with too few genes
go.env <- go.env[unlist(lapply(go.env, length))>5]
# append descriptions to the GO names
desc <- unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return("") } else { slot(x, "Term")}))
names(go.env) <- paste(names(go.env), desc)  # append description to the names
go.env <- list2env(go.env)  # convert to an environment

### Evaluate overdispersion of pre-defined gene sets
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = arg.n.cores, n.internal.shuffles = 0)
df.predefined <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
head(df.predefined)

### Evaluate overdispersion of 'de novo' gene sets
clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 150, n.cores = 1, plot = TRUE)
df.denovo <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
head(df.denovo)

dev.off()

### Visualize significant aspects of heterogeneity
# get full info on the top aspects
tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
# determine overall cell clustering
hc <- pagoda.cluster.cells(tam, varinfo)
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)

ccolors <- sampleTypeColor[as.character(myDesign[colnames(cd),"sampleType"])]
pdf(sprintf("%s.aspects.pdf",output.prefix),width = 20,height = 10)
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0, col.cols = ccolors, top=20)

col.cols <- rbind(groups = cutree(hc, 3))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = ccolors)

dev.off()

## compile a browsable app, showing top three clusters with the top color bar
#app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = col.cols, cell.clustering = hc, title = "NPCs")
# show app in the browser (port 1468)
#show.app(app, "SCS", browse = TRUE, port = 1468)  
#saveRDS(app, file = "SCS.app.rds")

save.image(file = output.rdata)

top.pathway <- sub("#PC1# ","",names(tamr2$cnam))
top.pathway.pattern <- sapply(top.pathway,function(p){ 
       pdf(sprintf("%s.pathway.%s.pdf",output.prefix,gsub("[:/]","_",p)),width = 20,height = 10)
       tryCatch(pagoda.show.pathways(p, varinfo, go.env, cell.clustering = hc, margins = c(1,5), n.genes = 30, show.cell.dendrogram = TRUE, colcols = ccolors, showRowLabels = TRUE, showPC = TRUE, main=p), 
                error = function(e) e, 
                finally = loginfo(sprintf("pathway finished: (%s)",p)) )
       dev.off()
        
})

# get cell cycle signature and view the top genes
#cc.pattern <- pagoda.show.pathways(c("GO:0000278 mitotic cell cycle","GO:0000280 nuclear division", "GO:0007067 mitotic nuclear division"), varinfo, go.env, show.cell.dendrogram = TRUE, cell.clustering = hc, showRowLabels = TRUE)
# subtract the pattern
#varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)
#pwpca.cc <- pagoda.pathway.wPCA(varinfo.cc, go.env, n.components = 1, n.cores = arg.n.cores, n.internal.shuffles = 0)
#df.predefined.cc <- pagoda.top.aspects(pwpca.cc, return.table = TRUE, plot = TRUE, z.score = 1.96)
#head(df.predefined.cc)

save.image(file = output.rdata)
