#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-e", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
###parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
parser$add_argument("-y", "--onlyY", action="store_true", default=FALSE, help="only center the data [default %(default)s]")
parser$add_argument("-k", "--notSC3", action="store_true", default=FALSE, help="don't run SC3 pipeline [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
geneFile <- args$geneFile
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
out.dir <- args$outDir
sample.id <- args$sample
args.log <- args$log
args.center <- args$center
nKeep <- args$nKeep
args.notSC3 <- args$notSC3
args.notFilter <- args$notFilter

##args <- commandArgs(T)
##if(length(args)<1){
##    cat("w.patch.P1116.R <nKeep>\n")
##    q()
##}
##nKeep=args[1]

#### TEST DATA 
###nKeep <- 1500
###designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/test/test.design.001.txt"
###inputFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/sf.byIndividual/sfDat.liver.txt.gz"
###cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
###clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/clonotype/clonotype.liver.txt"
###out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/test/OUT.centered"
###sample.id <- "liver"
###args.log <- TRUE
###args.center <- TRUE

dir.create(out.dir,recursive = T,showWarnings = F)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

##myDesign<-read.table(designFile,header=T,check.names=F,colClasses=c("factor","character","factor","factor"))
myDesign<-read.table(designFile,header=T,check.names=F,stringsAsFactors = F)
rownames(myDesign) <- myDesign$sample
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)

clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
patient.col.list <- patientColorListFromMyDesign(myDesign)

suppressPackageStartupMessages(library("R.utils"))
if(grepl("\\.scran\\.RData$",inputFile,perl=T)){
    lenv <- loadToEnv(inputFile)
    Y <- exprs(lenv[["sce.norm"]])
    args.notFilter <- T
    args.log <- F
    args.center <- F
}else if(grepl("RData$",inputFile,perl=T)){
    lenv <- loadToEnv(inputFile)
    Y <- lenv[["Y"]]
    ###obj.scdn <- lenv[["obj.scdn"]]
    ###Y <- obj.scdn@normalized.endo
}else{
    in.table <- read.table(inputFile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
}

sname <- intersect(rownames(myDesign),colnames(Y))
myDesign <- myDesign[sname,,drop=F]
Y <- Y[,sname,drop=F]

if(!args.notFilter){
    f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
    Y <- Y[f,]
}
if(!is.null(geneFile) && file.exists(geneFile)){
    geneTable <- read.table(geneFile,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
    Y <- Y[unique(as.character(geneTable$geneID)),]
}

if(args.log) { Y <- log2(Y+1) }
if(args.center){
    Y.new <- c()
    for(pp in unique(myDesign$patient)){
        Y.block <- t(scale(t(Y[,subset(myDesign,patient==pp,select="sample",drop=T)]),center = T,scale = F))
        Y.new <- cbind(Y.new,Y.block)
        print(apply(Y.block[1:4,],1,mean))
        print(apply(Y.block[1:4,],1,sd))
    }
    Y <- Y.new
    cat("test: all.equal(colnames(Y),rownames(myDesign))\n")
    print(all.equal(colnames(Y),rownames(myDesign)))
    Y <- Y[,sname,drop=F]
    cat("test: all.equal(colnames(Y),rownames(myDesign))\n")
    print(all.equal(colnames(Y),rownames(myDesign)))
}
#if(args$onlyY){
    save(Y,file=sprintf("%s/%s.Y.RData",out.dir,sample.id))
    loginfo("Y data saved.")
#    q()
#}

loginfo("... all samples.")

doit <- function(dat.plot,extra="",myseed=NULL){
    pca.res <- runPCAAnalysis(dat.plot,sprintf("%s/%s%s.het.PCA",out.dir,sample.id,extra),
                   myDesign[,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                   ntop=NULL,main=sample.id)

    tsne.res <- runTSNEAnalysis(dat.plot,sprintf("%s/%s%s.het.tSNE",out.dir,sample.id,extra),
                    col.points = sampleTypeColor[as.character(myDesign$sampleType)],
                    legend=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                             unique(myDesign$libType)),
                    col.legend=c(sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                                 rep("black",length(unique(myDesign$libType)))),
                    pch=(as.numeric(as.factor(myDesign$libType))-1) %% 26,cex=0.7,
                    pch.legend=c(rep(16,sum(names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType)))),
                                 (seq_along(unique(myDesign$libType))-1) %% 26),
                    myseed=myseed
                    )

    ##for(nn in unique(c(50,100,150,200,250,300,350,400,450,500,1000,nKeep,2000,100000)))
    for(nn in unique(c(50,100)))
    {
    runHierarchicalClusteringAnalysis(dat.plot,
                    sprintf("%s/%s%s.het.hclustering.n%s",out.dir,sample.id,extra,ifelse(nn>99999,"All",nn)),
                    myDesign[sname,"sampleType"],
                    sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                    clonotype.col=clonotype.strict.data,
                    patient.col.list = patient.col.list,
                    ntop=nn,
                    complexHeatmap.use=TRUE,do.cuttree=T,
                    verbose=TRUE,main="Variable Genes")
    }
    if(!args.notSC3){
        sc3.res <- runSC3Analysis(dat.plot,sprintf("%s/%s%s.SC3",out.dir,sample.id,extra),
                       myDesign[sname,"sampleType"],
                       sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
                       do.log.scale=FALSE,n.cores=10)
    }
    ###save(sc3.res,file=sprintf("%s/%s.SC3.RData",out.dir,sample.id))
}

Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=entrezToXXX(names(Y.sd.sort)))
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s/%s.Y.sd.sort.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)

g.f <- head(names(Y.sd.sort),n=nKeep)
doit(Y[g.f,],extra=".centered")
###doit(Y[names(Y.sd.sort)[1:nKeep],],extra=".centered.correct",myseed = 1471060113)
###doit(Y[names(Y.sd.sort)[1:nKeep],],extra=".centered.problem",myseed = 1471060113)

#save.image(sprintf("%s/%s.all.RData",out.dir,sample.id))

