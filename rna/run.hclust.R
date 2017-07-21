#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
#suppressPackageStartupMessages(require("gplots"))
#suppressPackageStartupMessages(require("ComplexHeatmap"))
#suppressPackageStartupMessages(require("circlize"))
#suppressPackageStartupMessages(require("gridBase"))
suppressPackageStartupMessages(require("dendextend"))
#suppressPackageStartupMessages(require("RColorBrewer"))


parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
#parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", type="character", required=TRUE, help="cellTypeColorFile")
parser$add_argument("-y", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-t", "--geneListFile", type="character", help="gene list file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-n", "--colname", type="character", default="sampleType", help="column name used to lable samples  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, 
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-k", "--nKeep", type="integer", default=3000, help="number of genes used [default %(default)s]")
parser$add_argument("-e", "--seed", type="integer", help="random number seed ")
args <- parser$parse_args()

print(args)

out.dir <- args$outDir
sample.id <- args$sample
designFile <- args$designFile
#cellTypeColorFile <- args$cellTypeColorFile
inputFile <- args$inputFile
args.log <- args$log
args.center <- args$center
args.colname <- args$colname
clonotype.file <- args$clonotypeFile
nKeep <- args$nKeep
rn.seed <- args$seed
geneListFile <- args$geneListFile

out.prefix <- sprintf("%s/%s.%s",out.dir,sample.id,"hclust")
#### TEST DATA P0205
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/test"
#sample.id <- "P1116"
#designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/OUT.scLVM.byTCellType/P1116/sfIgnoreERCC/TC/P1116.designUsed.txt"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#inputFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/OUT.scLVM.byTCellType/P1116/sfIgnoreERCC/TC/P1116.het.countGeneData.sfNormalized"
#args.log <- TRUE
#clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0729/filtered_TCR_summary/P0729.summary.cell.reassigneClonotype.txt"

dir.create(out.dir,recursive = T,showWarnings = F)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

##myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F)
##sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
sampleTypeColor <- auto.colSet(n=nlevels(myDesign[,args.colname]),name="Set1")
names(sampleTypeColor) <- levels(myDesign[,args.colname])

if(grepl("RData$",inputFile,perl=T)){
    suppressPackageStartupMessages(library("R.utils"))
    lenv <- loadToEnv(inputFile)
    Y <- lenv[["Y"]]
}else{
    in.table <- read.table(inputFile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
}
sname <- intersect(rownames(myDesign),colnames(Y))
myDesign <- myDesign[sname,,drop=F]
Y <- Y[,sname,drop=F]
loginfo(sprintf("before simple filter: %d genes",nrow(Y)))
f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
Y <- Y[f,]
loginfo(sprintf("after simple filter: %d genes",nrow(Y)))
if(args.log) { Y <- log2(Y+1) }
if(args.center){
    Y.new <- c()
    for(pp in unique(myDesign$patient)){
        Y.block <- t(scale(t(Y[,subset(myDesign,patient==pp,select="sample",drop=T)]),center = T,scale = F))
        Y.new <- cbind(Y.new,Y.block)
    }
    Y <- Y.new
    Y <- Y[,sname,drop=F]
}
print(dim(Y))

g.GNAME <- entrezToXXX(rownames(Y))
names(g.GNAME) <- rownames(Y)

clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
patient.col.list <- patientColorListFromMyDesign(myDesign)

Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=entrezToXXX(names(Y.sd.sort)))
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s.%s.Y.sd.sort.txt",out.prefix,sample.id),row.names = F,sep = "\t",quote = F)
if(!is.null(geneListFile) && file.exists(geneListFile)){
    geneList <- read.table(geneListFile,header = T,check.names = F,stringsAsFactors = F)
    g.f <- intersect(geneList$geneID,rownames(Y))
}else{
    g.f <- head(names(Y.sd.sort),n=nKeep)
}

dat.plot <- Y[g.f,]
clustering.distance <- "spearman"
clustering.method <- "complete"
obj.hclust.col <- NULL
if(clustering.distance=="spearman" || clustering.distance=="pearson"){
    tryCatch({
            obj.hclust.col <- hclust(as.dist(1-cor(dat.plot,method=clustering.distance)),method=clustering.method)
            ###branch.col <- color_branches(as.dendrogram(obj.hclust.col),k=k.col)
            },
            error = function(e) { 
                cat("using spearman or pearson as distance failed; will try to fall back to use euler distance ... \n"); e 
            }
        )
}
if(is.null(obj.hclust.col)){
    tryCatch({
        obj.hclust.col <- hclust(dist(t(dat.plot)),method=clustering.method)
    },
    error = function(e) {
        cat("failed using euler distance ... \n"); e 
    }
}
dend <- as.dendrogram(obj.hclust.col)
print(sampleTypeColor[myDesign[labels(dend),args.colname]])

pdf(sprintf("%s.pdf",out.prefix),width=12,height=8)
par(mar=c(10,8,4,2),cex.main=1.5,cex.lab=1.5)
plot(dend,main=sample.id,sub="",xlab="",cex=1.0*50/max(ncol(dat.plot),32))
colored_bars(colors = sampleTypeColor[myDesign[labels(dend),args.colname]], dend = dend, 
             rowLabels=args.colname, cex.rowLabels=1.1, y_scale=0.1)
legend("topright",legend=names(sampleTypeColor),fill=sampleTypeColor,border=sampleTypeColor,cex=1.5)
dev.off()

