#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-z", "--clusterColorFile", type="character", 
                    help="clusterColorFile [default %(default)s]")
##parser$add_argument("-y", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-t", "--geneListFile", type="character", help="gene list file")
parser$add_argument("-j", "--geneOnPCA", type="character", help="gene list file contain genes will be showed on PCA")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-q", "--save", action="store_true", default=FALSE, help="save object in .RData file [default %(default)s]")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-n", "--colname", type="character", default="sampleType", help="column name used to lable samples  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-a", "--rename", action="store_true", default=FALSE, help="rename sampleType [default %(default)s]")
parser$add_argument("-b", "--dbscan", action="store_true", default=FALSE, help="dbscan [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, 
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-k", "--nKeep", type="integer", default=3000, help="number of genes used [default %(default)s]")
parser$add_argument("-e", "--seed", type="integer", help="random number seed ")
parser$add_argument("-u", "--measure", type="character", help="measure (tpm)  [default %(default)s]")
parser$add_argument("-w", "--kbatch", type="character", default="2,3,4,5,6,7,8,9,10", help="kbatch [default %(default)s]")
args <- parser$parse_args()

print(args)

out.dir <- args$outDir
sample.id <- args$sample
designFile <- args$designFile
cellTypeColorFile <- args$cellTypeColorFile
clusterColorFile <- args$clusterColorFile
inputFile <- args$inputFile
args.log <- args$log
args.center <- args$center
args.colname <- args$colname
###clonotype.file <- args$clonotypeFile
nKeep <- args$nKeep
rn.seed <- args$seed
geneListFile <- args$geneListFile
args.save <- args$save
args.rename <- args$rename
args.dbscan <- args$dbscan
args.geneOnPCA <- args$geneOnPCA
args.measure <- args$measure
args.kbatch <- args$kbatch
if(!is.null(args.kbatch)) { args.kbatch <- as.numeric(unlist(strsplit(args.kbatch,","))) }

out.prefix <- sprintf("%s/%s.%s",out.dir,sample.id,"tSNEandPCA")

dir.create(out.dir,recursive = T,showWarnings = F)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("gplots"))
suppressPackageStartupMessages(require("factoextra"))
suppressPackageStartupMessages(require("FactoMineR"))
suppressPackageStartupMessages(require("fields"))

#####
g.inputD <- processInput(designFile,cellTypeColorFile = cellTypeColorFile,inputFile,args.notFilter=F,geneFile = NULL,
                         args.center,args.log,args.norm.exprs = F,args.measure = args.measure)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
##### 

### rename "sampleType"
if(args.rename){
    myDesign$sampleType <- paste(myDesign$stype,sapply(strsplit(myDesign$sampleType,""),function(x){x[1]}),sep=".")
}
## check majorCluster
if(!"majorCluster" %in% colnames(myDesign)){
    myDesign[,"majorCluster"] <- myDesign[,"sampleType"]
}

#### myajorCluster color
nMajor <- length(unique(myDesign$majorCluster))
if(!is.null(clusterColorFile) && file.exists(clusterColorFile)){
    majorClusterColor <- read.SampleTypeColor(clusterColorFile)
    majorClusterColor <- majorClusterColor[ names(majorClusterColor) %in% unique(myDesign[colnames(Y),"majorCluster"]) ]
}else{
    majorClusterColor <- structure(colorRampPalette(brewer.pal(nMajor,"Paired"))(nMajor),
                              names=unique(myDesign$majorCluster))
    majorClusterColor.df <- data.frame(sampleType=names(majorClusterColor),color=majorClusterColor)
    write.table(majorClusterColor.df,sprintf("%s.%s.majorClusterColor.txt",out.prefix,sample.id),row.names = F,sep = "\t",quote = F)
}
print("sampleTypeColor")
print(sampleTypeColor)
print("majorClusterColor")
print(majorClusterColor)

print(dim(Y))
print(Y[1:4,1:6])

if(!is.null(geneListFile) && file.exists(geneListFile)){
    ###geneList <- read.table(geneListFile,sep="\t",header = T,check.names = F,stringsAsFactors = F)
    geneList <- read.delim(geneListFile,sep="\t",header = T,check.names = F,stringsAsFactors = F)
    geneList <- head(geneList,n=nKeep)
    g.f <- intersect(as.character(geneList$geneID),rownames(Y))
}else{
    Y.sd <- apply(Y, 1, sd)
    Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
    #Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=entrezToXXX(names(Y.sd.sort)))
    Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=g.GNAME[names(Y.sd.sort)])
    Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
    write.table(Y.sd.df,sprintf("%s.%s.Y.sd.sort.txt",out.prefix,sample.id),row.names = F,sep = "\t",quote = F)
    g.f <- head(names(Y.sd.sort),n=nKeep)
}

loginfo("tSNE begin.")
tsne.out <- NULL
tsne.dbscan.out <- NULL
tsne.out <- runTSNEAnalysis(Y[g.f,],out.prefix,
                            legend=names(majorClusterColor),
                            col.points=majorClusterColor[as.character(myDesign$majorCluster)],
                            col.legend=majorClusterColor,
                            pch=20,pch.legend=20,
                            inPDF=TRUE,eps.clus = args.kbatch,dims=2,k=NULL,
                            do.dbscan=F,myseed=rn.seed)

tsne.dbscan.out <- runTSNEAnalysis(Y[g.f,],sprintf("%s.dbScan",out.prefix),
                                   legend=names(sampleTypeColor),
                                   col.points=sampleTypeColor[myDesign$sampleType],
                                   col.legend=sampleTypeColor,
                                   pch=20,pch.legend=20,
                                   inPDF=TRUE,eps=args.kbatch,dims=2,k=NULL,
                                   do.dbscan=args.dbscan,myseed=rn.seed,preSNE=tsne.out)

                                   ###do.dbscan=T,myseed=rn.seed,preSNE=tsne.out[["30"]]$Rtsne.res)
if(args.save){ save(tsne.out,tsne.dbscan.out,myDesign,
                    sampleTypeColor,majorClusterColor,
                    file=sprintf("%s.obj.RData",out.prefix)) }
loginfo("tSNE done.")

loginfo("PCA begin.")
pca.res.sampleType <- runPCAAnalysis(Y[g.f,],sprintf("%s/%s.sampleType.PCA",out.dir,sample.id),
               myDesign[,"sampleType"],
               sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[,"sampleType"]))],
               ntop=NULL,main=sample.id,myseed=rn.seed)
pca.res.majorCluster <- runPCAAnalysis(Y[g.f,],sprintf("%s/%s.majorCluster.PCA",out.dir,sample.id),
               myDesign[,"majorCluster"],
               majorClusterColor[names(majorClusterColor) %in% unique(as.character(myDesign[,"majorCluster"]))],
               ntop=NULL,main=sample.id,myseed=rn.seed)
loginfo("PCA done.")

if(!is.null(args.geneOnPCA) && file.exists(args.geneOnPCA)){
    gene.to.show <- read.table(args.geneOnPCA,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    gene.to.show <- structure(as.character(gene.to.show$geneID),names=gene.to.show$symbol)
    dir.create(sprintf("%s.perGene.PCA",out.prefix),showWarnings = F,recursive = T)
    dir.create(sprintf("%s.perGene.tSNE",out.prefix),showWarnings = F,recursive = T)
    for(i in seq_along(gene.to.show)){
        gid <- gene.to.show[i]
        gname <- names(gene.to.show)[i]
        s.f <- order(Y[gid,],decreasing=F)
        Y.level <- pretty(Y[gid,],n=10)
        ppalette <- brewer.pal(9,"YlOrRd")
        gid.color <- as.character(cut(Y[gid,s.f],
                                      breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
                                      labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
        ####### PCA #######
        pdf(sprintf("%s.perGene.PCA/PCA.%s.pdf",out.prefix,gname),width = 9,height = 8)
        par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=3.0,cex.axis=1.5)
        plot(pca.res.sampleType$ind$coord[,"Dim.1"],pca.res.sampleType$ind$coord[,"Dim.2"],
             t='n',pch=16,col="lightgray", main=sprintf("%s",gname),xlab="Dim1",ylab="Dim2")
        points(pca.res.sampleType$ind$coord[s.f,"Dim.1"],pca.res.sampleType$ind$coord[s.f,"Dim.2"],col=gid.color,pch=16)
        image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),legend.only=TRUE,col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=5.0)
        dev.off()
        ####### tSNE #######
        slist.tsneUsed <- rownames(tsne.out$`30`$Rtsne.res$Y)
        i.Y <- Y[gid,slist.tsneUsed]
        s.f <- order(i.Y,decreasing=F)
        gid.color <- as.character(cut(i.Y[s.f],
                                      breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
                                      labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
        pdf(sprintf("%s.perGene.tSNE/tSNE.%s.pdf",out.prefix,gname),width = 9,height = 8)
        par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=3.0,cex.axis=1.5)
        plot(tsne.out$`30`$Rtsne.res$Y,
             t='n',pch=16,col="lightgray", main=sprintf("%s",gname),xlab="Dim1",ylab="Dim2")
        points(tsne.out$`30`$Rtsne.res$Y[s.f,],col=gid.color,pch=16)
        image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),legend.only=TRUE,col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=5.0)
        dev.off()
    }

}


