#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile")
parser$add_argument("-e", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
args <- parser$parse_args()

print(args)

out.dir <- args$outDir
sample.id <- args$sample
designFile <- args$designFile
cellTypeColorFile <- args$cellTypeColorFile
inputFile <- args$inputFile
args.log <- args$log
clonotype.file <- args$clonotypeFile

#### TEST DATA 
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/test"
#sample.id <- "P1116"
#designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/OUT.scLVM.byTCellType/P1116/sfIgnoreERCC/TC/P1116.designUsed.txt"
#designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/OUT.scLVM.byTCellType/P1116/sfEndo/TC/P1116.designUsed.txt"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#inputFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/OUT.scLVM.byTCellType/P1116/sfIgnoreERCC/TC/P1116.het.countGeneData.sfNormalized"
#inputFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/OUT.scLVM.byTCellType/P1116/sfEndo/TC/P1116.RData"
#args.log <- TRUE
#clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P1116/filtered_TCR_summary/P1116.summary.cell.reassigneClonotype.txt"

dir.create(out.dir,recursive = T,showWarnings = F)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)

clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
patient.col.list <- patientColorListFromMyDesign(myDesign)

if(grepl("RData$",inputFile,perl=T)){
    suppressPackageStartupMessages(library("R.utils"))
    lenv <- loadToEnv(inputFile)
    obj.scdn <- lenv[["obj.scdn"]]
    Y <- rbind(obj.scdn@normalized.endo,obj.scdn@normalized.ERCC)
}else{
    in.table <- read.table(inputFile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
}
f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
Y <- Y[f,]
cIdx <- grep(pattern = "^ERCC",rownames(Y),perl=T)
if(args.log) { Y <- log2(Y+1) }

loginfo("... all samples.")

sname <- intersect(rownames(myDesign),colnames(Y))
myDesign <- myDesign[sname,,drop=F]
Y <- Y[,sname,drop=F]

run.RUVnormalize <- function(Y,cIdx){
    require(RUVnormalize)
    require(spams)
    ret.list <- list()
    m <- nrow(Y)
    n <- ncol(Y)
    #p <- ncol(X)
    p <- 1

    ## Naive RUV-2 no shrinkage
    k <- 20
    nu <- 0
    nsY <- naiveRandRUV(Y, cIdx, nu.coeff=0, k=k)
    ret.list[["nsY"]] <- nsY
    ## Naive RUV-2 + shrinkage
    k <- m
    nu.coeff <- 1e-3
    nY <- naiveRandRUV(Y, cIdx, nu.coeff=nu.coeff, k=k)
    ret.list[["nY"]] <- nY
    ## Iterated ridge
    cEps <- 1e-6
    maxIter <- 30
    paramXb <- list()
    paramXb$K <- p
    paramXb$D <- matrix(c(0.),nrow = 0,ncol=0)
    paramXb$batch <- TRUE
    paramXb$iter <- 1
    paramXb$mode <- 'PENALTY' #2
    paramXb$lambda <- 6e-2
    paramXb$lambda2 <- 0
    iRes <- iterativeRUV(Y, cIdx, scIdx=NULL, paramXb, k=nrow(Y), nu.coeff=1e-3/2, cEps, maxIter, Wmethod='svd', wUpdate=11)
    nrcY <- iRes$cY
    ret.list[["nrcY"]] <- nrcY

    return(ret.list)
}
#RUVnormalize.out <- run.RUVnormalize(t(Y[,1:100]),cIdx)
RUVnormalize.out <- run.RUVnormalize(t(Y),cIdx)
my.nY <- t(RUVnormalize.out$nY)
my.nY.sd <- apply(my.nY, 1, sd)
my.nY.sd.sidx <- sort(my.nY.sd,decreasing=TRUE,index.return=TRUE)$ix
nKeep <- 1500

### output normalized result
out.Y.df <- data.frame(geneID=rownames(my.nY),geneSymbol=entrezToXXX(rownames(my.nY)))
out.Y.df <- cbind(out.Y.df,my.nY)
f <- grepl("^ERCC",out.Y.df$geneID,perl=T)
out.Y.df <- out.Y.df[!f,]
write.table(out.Y.df,file = sprintf("%s/%s.RUVnormalized.nY.txt",out.dir,sample.id), 
            quote = F,sep = "\t",row.names = F,col.names = T)

###
doit <- function(dat.plot,extra=""){
    pca.res <- runPCAAnalysis(dat.plot,sprintf("%s/%s%s.het.PCA",out.dir,sample.id,extra),
                   myDesign[,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                   ntop=NULL,main=sample.id)

    tsne.res <- runTSNEAnalysis(dat.plot,sprintf("%s/%s%s.het.tSNE",out.dir,sample.id,extra),
                    col.points = sampleTypeColor[as.character(myDesign$sampleType)],
                    legend=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],levels(myDesign$libType)),
                    col.legend=c(sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],rep("black",length(levels(myDesign$libType)))),
                    pch=(as.numeric(myDesign$libType)-1) %% 26,cex=0.7,
                    pch.legend=c(rep(16,sum(names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType)))),(seq_along(levels(myDesign$libType))-1) %% 26)
                    )

    for(nn in c(50,100,150,200,250,300,350,400,450,500,1000,2000,100000))
    {
    runHierarchicalClusteringAnalysis(dat.plot,
                    sprintf("%s/%s%s.het.hclustering.n%s",out.dir,sample.id,extra,ifelse(nn>99999,"All",nn)),
                    myDesign[sname,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                    clonotype.col=clonotype.strict.data,
                    patient.col.list = patient.col.list,
                    ntop=nn,
                    complexHeatmap.use=TRUE,do.cuttree=T,
                    verbose=FALSE,main="Variable Genes")
    }

    sc3.res <- runSC3Analysis(dat.plot,sprintf("%s/%s%s.SC3",out.dir,sample.id,extra),
                   myDesign[sname,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
                   do.log.scale=FALSE,n.cores=10)
    ###save(sc3.res,file=sprintf("%s/%s.SC3.RData",out.dir,sample.id))
}

doit(my.nY[my.nY.sd.sidx[1:nKeep],],extra="")

Y.sd <- apply(Y, 1, sd)
Y.sd.sidx <- sort(Y.sd,decreasing=TRUE,index.return=TRUE)$ix
doit(Y[Y.sd.sidx[1:nKeep],],extra=".uncorrected")

my.nrcY <- t(RUVnormalize.out$nrcY)
my.nrcY.sd <- apply(my.nrcY, 1, sd)
my.nrcY.sd.sidx <- sort(my.nrcY.sd,decreasing=TRUE,index.return=TRUE)$ix
doit(my.nrcY[my.nrcY.sd.sidx[1:nKeep],],extra=".RUVIter")

save.image(sprintf("%s/%s.all.RData",out.dir,sample.id))

