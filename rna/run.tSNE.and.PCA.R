#!/usr/bin/env Rscript

uppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", type="character", required=TRUE, help="cellTypeColorFile")
parser$add_argument("-e", "--clonotypeFile", type="character", help="clonotype file")
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
cellTypeColorFile <- args$cellTypeColorFile
inputFile <- args$inputFile
args.log <- args$log
args.center <- args$center
args.colname <- args$colname
clonotype.file <- args$clonotypeFile
nKeep <- args$nKeep
rn.seed <- args$seed

out.prefix <- sprintf("%s/%s.%s",out.dir,sample.id,"tSNEandPCA")
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
print(Y[1:4,1:6])

g.GNAME <- entrezToXXX(rownames(Y))
names(g.GNAME) <- rownames(Y)

clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
patient.col.list <- patientColorListFromMyDesign(myDesign)

Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=entrezToXXX(names(Y.sd.sort)))
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s.%s.Y.sd.sort.txt",out.prefix,sample.id),row.names = F,sep = "\t",quote = F)

loginfo("tSNE begin.")
tsne.out <- runTSNEAnalysis(Y[head(names(Y.sd.sort),n=nKeep),],out.prefix,
                            legend=names(sampleTypeColor),
                            col.points=sampleTypeColor[myDesign[,args.colname]],
                            col.legend=sampleTypeColor,
                            pch=16,pch.legend=16,
                            inPDF=TRUE,eps=2.0,dims=2,k=NULL,
                            do.dbscan=F,myseed=rn.seed,width.pdf = 14,margin.r = 24,legend.inset = -0.45)

#tsne.dbscan.out <- runTSNEAnalysis(Y[head(names(Y.sd.sort),n=nKeep),],sprintf("%s.dbScan",out.prefix),
#                                   legend=names(sampleTypeColor),
#                                   col.points=sampleTypeColor[sample.desc$sampleType],
#                                   col.legend=sampleTypeColor,
#                                   pch=16,pch.legend=16,
#                                   inPDF=TRUE,eps=1.00,dims=2,k=NULL,
#                                   do.dbscan=T,myseed=rn.seed,preSNE=tsne.out[["30"]]$Rtsne.res,width.pdf = 10,
#                                   margin.r = 8,legend.inset = -0.21)

if(args.save){ save(tsne.out,tsne.dbscan.out,sample.desc,file=sprintf("%s.obj.RData",out.prefix)) }
loginfo("tSNE done.")

pca.res <- runPCAAnalysis(Y[head(names(Y.sd.sort),n=nKeep),],sprintf("%s/%s.het.PCA",out.dir,sample.id),
               myDesign[,args.colname],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[,args.colname]))],
               ntop=NULL,main=sample.id)

for(nn in c(50,100,150,200,250,300,350,400,450,500,1000,2000,3000,100000))
{
    runHierarchicalClusteringAnalysis(Y,
                    sprintf("%s/%s.het.hclustering.n%s",out.dir,sample.id,ifelse(nn>99999,"All",nn)),
                    myDesign[,args.colname],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                    clonotype.col=clonotype.strict.data,
                    patient.col.list = patient.col.list,
                    ntop=nn,
                    complexHeatmap.use=TRUE,do.cuttree=T,
                    verbose=TRUE,main="Variable Genes")
}

