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
parser$add_argument("-a", "--rename", action="store_true", default=FALSE, help="rename sampleType (need stype field in designFile) [default %(default)s]")
parser$add_argument("-b", "--dbscan", action="store_true", default=FALSE, help="dbscan [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-f", "--noFilter", action="store_true", default=FALSE, help="don't apply filter on genes [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, 
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-x", "--scale", action="store_true", default=FALSE, help="zscore in geneOnTSNE plot [default %(default)s]")
parser$add_argument("-k", "--nKeep", type="integer", default=3000, help="number of genes used [default %(default)s]")
parser$add_argument("-e", "--seed", type="integer", help="random number seed ")
parser$add_argument("-u", "--measure", type="character", help="measure (exprs ,tpm, norm_exprs)  [default %(default)s]")
parser$add_argument("-m", "--mark", type="character", help="mark points by columns in design file [default %(default)s]")
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
args.mark <- args$mark
if(!is.null(args.mark)) { args.mark <- unlist(strsplit(args.mark,",")) }
args.scale <- args$scale
args.noFilter <- args$noFilter

out.prefix <- sprintf("%s/%s.%s",out.dir,sample.id,"tSNEandPCA")

dir.create(out.dir,recursive = T,showWarnings = F)
source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("factoextra"))
suppressPackageStartupMessages(library("FactoMineR"))
suppressPackageStartupMessages(library("fields"))

mark.plot <- function(x,mark.name,d,out.file,x.dim=c(1,2),colSet=NULL){
    ss <- intersect(rownames(x),rownames(d))
    x <- x[ss,,drop=F]
    d <- d[ss,,drop=F]
    if(class(d[,mark.name])=="character" || class(d[,mark.name])=="factor" || class(d[,mark.name])=="logical"){
        nn <- sort(unique(as.character((d[,mark.name]))))
        f.nn <- is.na(nn) | nn=="NA" | nn=="."
        if(is.null(colSet)){
            colSet <- structure(auto.colSet(n = length(nn[!f.nn]),name = "Dark2"),
                                names=nn[!f.nn])
            if(sum(f.nn)>0){
                colSet <- c(colSet,structure(rep("lightgray",sum(f.nn)),names=nn[f.nn]))
            }
        }
        print("XXXX")
        print(nn)
        print(colSet)
        plot.tsne.points(x,out.file,tsne.col.points=colSet[d[rownames(x),mark.name]],
                         col.tsne.legend=colSet,
                         tsne.legend=names(colSet),pch=20,nclusters=length(nn),peak=NULL,main=mark.name,x.dim=x.dim)
    }else if(class(d[,mark.name])=="integer" || class(d[,mark.name])=="numeric"){
        x <- x[,x.dim,drop=F]
        s.f <- order(d[,mark.name],decreasing=F)
        Y.level <- pretty(d[,mark.name],n=10)
        ppalette <- brewer.pal(9,"YlGnBu")
        Y.color <- as.character(cut(d[s.f,mark.name],
                                      breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
                                      labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
        pdf(out.file,width=9,height=8)
        par(mar=c(5,5,4,9),cex.lab=1.8,cex.main=2.5,cex.axis=1.5)
        plot(x[s.f,1],x[s.f,2],
             t='n',pch=20,col="lightgray", main=sprintf("%s",mark.name),
             xlab=sprintf("Dim%d",x.dim[1]),
             ylab=sprintf("Dim%d",x.dim[2]))
        points(x[s.f,1],x[s.f,2],col=Y.color,pch=20,cex=auto.point.size(nrow(x)))
        image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),
                   legend.only=TRUE,col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=8.0)
        dev.off()
    }
}

#####
g.inputD <- processInput(designFile,cellTypeColorFile = cellTypeColorFile,inputFile,args.notFilter=args.noFilter,geneFile = NULL,
                         args.center,args.log,args.norm.exprs = F,args.measure = args.measure)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
##### 

#print("sampleTypeColor (00)")
#print(sampleTypeColor)
#print(unique(as.character(myDesign[,"sampleType"])))
### rename "sampleType"
if(args.rename){
    myDesign$sampleType <- paste(myDesign$stype,sapply(strsplit(myDesign$sampleType,""),function(x){x[1]}),sep=".")
    #### reset sampleTypeColor
    if(!is.null(cellTypeColorFile) && file.exists(cellTypeColorFile)){
        sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
        sampleTypeColor <- sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[,"sampleType"]))]
    }else{
        sampleTypeColor <- auto.colSet(n=length(unique(myDesign[,"sampleType"])),name="Paired")
        names(sampleTypeColor) <- unique(myDesign[,"sampleType"])
    }
}

## check majorCluster
majorClusterColor <- NULL
if("majorCluster" %in% colnames(myDesign)){
    #### myajorCluster color
    nMajor <- length(unique(myDesign$majorCluster))
    if(!is.null(clusterColorFile) && file.exists(clusterColorFile)){
        majorClusterColor <- read.SampleTypeColor(clusterColorFile)
        majorClusterColor <- majorClusterColor[ names(majorClusterColor) %in% unique(myDesign[colnames(Y),"majorCluster"]) ]
    }else{
        majorClusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(nMajor),
                                  names=sort(unique(myDesign$majorCluster)))
        majorClusterColor.df <- data.frame(sampleType=names(majorClusterColor),color=majorClusterColor)
        write.table(majorClusterColor.df,sprintf("%s.%s.majorClusterColor.txt",out.prefix,sample.id),row.names = F,sep = "\t",quote = F)
    }
    args.mark <- unique(c(args.mark,"majorCluster"))
}
print("sampleTypeColor")
print(sampleTypeColor)
print("majorClusterColor")
print(majorClusterColor)

#print(dim(Y))
#print(Y[1:4,1:6])

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

#if(!is.null(args.preResult) && file.exists(args.preResult)){
#    penv <- loadToEnv(args.preResult)
#    tsne.out <- penv[["tsne.out"]]
#    tsne.dbscan.out <- penv[["tsne.dbscan.out"]]


tsne.out <- runTSNEAnalysis(Y[g.f,],out.prefix,
                            legend=names(sampleTypeColor),
                            col.points=sampleTypeColor[as.character(myDesign$sampleType)],
                            col.legend=sampleTypeColor,
                            pch=20,pch.legend=20,
                            inPDF=TRUE,eps.clus = args.kbatch,dims=2,k=NULL,
                            do.dbscan=F,myseed=rn.seed,do.scale=T)

for(mm in args.mark){
    mark.plot(tsne.out[["30"]]$Rtsne.res$Y,
              mark.name = mm,
              d = myDesign,
              out.file = sprintf("%s.TSNE.mark.%s.pdf",out.prefix,mm),
              colSet=if(mm=="majorCluster") majorClusterColor else NULL)
}

loginfo("tSNE done.")

loginfo("PCA begin.")
pca.res.sampleType <- runPCAAnalysis(Y[g.f,],sprintf("%s/%s.sampleType.PCA",out.dir,sample.id),
               myDesign[,"sampleType"],
               sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[,"sampleType"]))],
               ntop=NULL,main=sample.id,myseed=rn.seed,gid.mapping=g.GNAME)
##pca.res.majorCluster <- runPCAAnalysis(Y[g.f,],sprintf("%s/%s.majorCluster.PCA",out.dir,sample.id),
##               myDesign[,"majorCluster"],
##               majorClusterColor[names(majorClusterColor) %in% unique(as.character(myDesign[,"majorCluster"]))],
##               ntop=NULL,main=sample.id,myseed=rn.seed,gid.mapping=g.GNAME)

for(mm in args.mark){
    mark.plot(pca.res.sampleType$ind$coord,
              mark.name = mm,
              d = myDesign,
              out.file = sprintf("%s.PCA.mark.%s.PC1.PC2.pdf",out.prefix,mm),x.dim=c(1,2),
              colSet=if(mm=="majorCluster") majorClusterColor else NULL)
    mark.plot(pca.res.sampleType$ind$coord,
              mark.name = mm,
              d = myDesign,
              out.file = sprintf("%s.PCA.mark.%s.PC1.PC3.pdf",out.prefix,mm),x.dim=c(1,3),
              colSet=if(mm=="majorCluster") majorClusterColor else NULL)
    mark.plot(pca.res.sampleType$ind$coord,
              mark.name = mm,
              d = myDesign,
              out.file = sprintf("%s.PCA.mark.%s.PC2.PC3.pdf",out.prefix,mm),x.dim=c(2,3),
              colSet=if(mm=="majorCluster") majorClusterColor else NULL)
}

loginfo("PCA done.")
if(args.save){ save(tsne.out,myDesign,pca.res.sampleType,
                    sampleTypeColor,majorClusterColor,
                    file=sprintf("%s.obj.RData",out.prefix)) }

if(!is.null(args.geneOnPCA) && file.exists(args.geneOnPCA)){
    gene.to.show <- read.table(args.geneOnPCA,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    gene.to.show <- structure(as.character(gene.to.show[,1]),names=gene.to.show[,2])
    dir.create(sprintf("%s.perGene.PCA",out.prefix),showWarnings = F,recursive = T)
    dir.create(sprintf("%s.perGene.tSNE",out.prefix),showWarnings = F,recursive = T)
    for(i in seq_along(gene.to.show)){
        gid <- gene.to.show[i]
        gname <- names(gene.to.show)[i]
        if(!(gid %in% rownames(Y))){ 
            cat(sprintf("Gene not in data: %s (%s)\n",gid, gname))
            next
        }
        i.Y <- Y[gid,]
        s.f <- order(i.Y,decreasing=F)
        Y.level <- pretty(Y[gid,],n=10)
        ppalette <- brewer.pal(9,"YlOrRd")
        if(args.scale){
            i.Y <- scale(i.Y)
            #Y.level <- pretty(i.Y,n=10)
            Y.level <- c(-2.5,2.5)
            ppalette <- rev(brewer.pal(9,"RdBu"))
            print(gname)
            print(sprintf("mean: %4.4f, sd: %4.4f", mean(i.Y), sd(i.Y)))
            print(summary(i.Y))
            i.Y[i.Y>=2.5] <- 2.5
            i.Y[i.Y<=-2.5] <- -2.5
        }
        gid.color <- as.character(cut(i.Y[s.f],
                                      breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
                                      labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
        ####### PCA #######
        pdf(sprintf("%s.perGene.PCA/PCA.%s.pdf",out.prefix,gname),width = 9,height = 8)
        par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=3.0,cex.axis=1.5)
        plot(pca.res.sampleType$ind$coord[,"Dim.1"],pca.res.sampleType$ind$coord[,"Dim.2"],
             t='n',pch=16,col="lightgray", main=sprintf("%s",gname),xlab="Dim1",ylab="Dim2")
        points(pca.res.sampleType$ind$coord[s.f,"Dim.1"],pca.res.sampleType$ind$coord[s.f,"Dim.2"],col=gid.color,pch=16)
        image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),legend.only=TRUE,
                   col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=5.0)
        dev.off()
        ####### tSNE #######
        slist.tsneUsed <- rownames(tsne.out$`30`$Rtsne.res$Y)
        i.Y <- Y[gid,slist.tsneUsed]
        s.f <- order(i.Y,decreasing=F)
        if(args.scale){
            i.Y <- scale(i.Y)
            print(gname)
            print(sprintf("mean: %4.4f, sd: %4.4f", mean(i.Y), sd(i.Y)))
            print(summary(i.Y))
            i.Y[i.Y>=2.5] <- 2.5
            i.Y[i.Y<=-2.5] <- -2.5
            #Y.level <- pretty(i.Y,n=10)
            Y.level <- c(-2.5,2.5)
            ppalette <- rev(brewer.pal(9,"RdBu"))
        }
        gid.color <- as.character(cut(i.Y[s.f],
                                      breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
                                      labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
        
        pdf(sprintf("%s.perGene.tSNE/tSNE.%s.pdf",out.prefix,gname),width = 9,height = 8)
        par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=3.0,cex.axis=1.5)
        plot(tsne.out$`30`$Rtsne.res$Y,
             t='n',pch=16,col="lightgray", main=sprintf("%s",gname),xlab="Dim1",ylab="Dim2")
        points(tsne.out$`30`$Rtsne.res$Y[s.f,],col=gid.color,pch=16)
        image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),legend.only=TRUE,
                   col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=5.0)
        dev.off()
    }

}


