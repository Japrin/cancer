#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--in",required=T)
parser$add_argument("-o", "--out",required=T)
parser$add_argument("-a", "--a",default=0.05,type="double",help="downsample range's start point [default %(default)s]")
parser$add_argument("-b", "--b",default=0.05,type="double",help="downsample range's step [default %(default)s]")
parser$add_argument("-c", "--c",default=FALSE,action="store_true",
                    help="whether donwsample to predefined numbers [default %(default)s]")
parser$add_argument("-p", "--pair", action="store_true", default=TRUE, 
                    help="if specified, PE reads other wize SE reads [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose output [default %(default)s]")
parser$add_argument("-s", "--sample", default="SAMPLE", help="sample id [default %(default)s]")
parser$add_argument("-n", "--ncores", type="integer", default=4, 
                    help="n cores [default %(default)s]")
parser$add_argument("-m", "--multiple", type="integer", default=1, 
                    help="multiple sample [default %(default)s]")
parser$add_argument("-y", "--ylab",type="character",default="Number Of Detected Genes", help="ylab [default %(default)s]")
parser$add_argument("-f", "--genelist",type="character", help="gene list file, format: geneID[tab]geneSymbol [default %(default)s]")
parser$add_argument("-g", "--gfeature", default="/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/hg19.knownGene.RData", 
                    help="RData contains gene feature [default %(default)s]")
args <- parser$parse_args()

in.bam.file <- args[["in"]]
out.prefix <- args$out
paired_ends <- args$pair
verbose <- args$verbose
sample.id <- args$sample
gfeature.file <- args$gfeature
ds.start <- args$a
ds.step <- args$b
args.m <- args$multiple
n.cores <- args$ncores
args.c <- args$c
gene.list.file <- args$genelist
args.ylab <- args$ylab

#### For test only
##in.bam.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/OUT/P0205/TTS232-0205/bams/TTS232-0205.analyzed.bam"
##out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/liver/preliminary/test.geneNumber_saturation"
##paired_ends <- TRUE
##verbose <- TRUE
##sample.id <- "SAMPLE"
##gfeature.file <- "/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/hg19.knownGene.RData"
##ds.start <- 0.2
##ds.step <- 0.2
##args.m <- 10
##n.cores <- 4

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

#### functin definition
suppressPackageStartupMessages(library("HTSeqGenie"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(require("plyr"))
suppressPackageStartupMessages(require("doParallel"))
suppressPackageStartupMessages(require("plotrix"))
suppressPackageStartupMessages(library("tictoc"))

computeWidth <- function(gf) {
  if (class(gf)=="GRangesList") sapply(width(reduce(gf)), sum)
  else width(gf)
}

rpkm <- function (counts, widths, nbreads)  {
  (1e9*counts) / (nbreads * as.numeric(widths))
}

tpm <- function (counts, widths)  {
  v <- counts / as.numeric(widths)
  s <- sum(v)
  (1e6*v) / s
}

countGenomicFeaturesSimple <- function(features,reads,nbanalyzedreads)
{
	overlaps <- suppressWarnings(findOverlaps(features, reads))
	nbreads.byfeature <- tabulate(queryHits(overlaps), nbins=queryLength(overlaps))
	nbreads <- length(unique(subjectHits(overlaps)))
	names(nbreads.byfeature ) <- names(features)

	widths <- computeWidth(features)
	rpkms <- rpkm(counts=nbreads.byfeature, widths=widths, nbanalyzedreads)
	tpms <- tpm(counts=nbreads.byfeature, widths=widths)
	df <- data.frame(name=names(nbreads.byfeature), count=nbreads.byfeature, width=widths, rpkm=rpkms, tpm=tpms)
	df
}


if(file.exists(sprintf("%s.RData",out.prefix))){
    lname <- load(sprintf("%s.RData",out.prefix))
}else{
    #### process input data
    if(!is.null(gene.list.file) && file.exists(gene.list.file)){
        gene.list <- read.table(gene.list.file,header = T,check.names = F,stringsAsFactors = F)
        gene.list$geneID <- as.character(gene.list$geneID)
        cat(sprintf("Total number genes in this list: %d (%s)\n",nrow(gene.list),gene.list.file))
    }else{
        gene.list <- NULL
    }

    reads <- readRNASeqEnds(in.bam.file,paired_ends=paired_ends)
    #head(reads)
    genomic_features <- get(load(gfeature.file))
    gc()
    nbanalyzedreads <- length(reads)
    detected.gn <- list()
    registerDoParallel(cores = n.cores)

    cat(sprintf("[%s] begin downsample & count reads:\n", Sys.time()))
    detected.gn <- llply(seq_len(args.m),function(ii){
        ds <- sample(1:nbanalyzedreads)
        ds.pct <- seq(ds.start,1,ds.step)
        if(args.c){
            ds.n <- c(100000,150000,200000,250000,seq(300000,3000000,200000))
            ds.n <- ds.n[ds.n <= nbanalyzedreads]
        }else{
            ds.n <- as.integer(nbanalyzedreads*ds.pct)
        }
        ds.res.byRPKM <- data.frame()
        ds.res.byTPM <- data.frame()
        ds.res.byCount <- data.frame() 
        tic(sprintf("one run (%d)",ii))
        hh <- sapply(seq_along(ds.n),function(i,n,pct){
            #tic(sprintf("countGenomicFeaturesSimple(%d reads)",n[i]))
            res.ds <- countGenomicFeaturesSimple(features = genomic_features$gene_exonic,
                                                 reads = reads[ds[1:n[i]]],nbanalyzedreads = n[i])
            #toc()
            if(ncol(ds.res.byRPKM)==0)
            {
                ds.res.byRPKM <<- res.ds[,c("name","width")]
                ds.res.byTPM <<- res.ds[,c("name","width")]
                ds.res.byCount <<- res.ds[,c("name","width")]
            }
            if(args.c){
                this.colname <- sprintf("%d",n[i])
            }else{
                this.colname <- paste0("",pct[i],"")
            }
            ds.res.byRPKM[,this.colname] <<- res.ds$rpkm
            ds.res.byTPM[,this.colname] <<- res.ds$tpm
            ds.res.byCount[,this.colname] <<- res.ds$count
            
        },ds.n,ds.pct)
        toc()
        e.TH <- c(0,0.1,0.5,1,2,3,5,10,30,50,100)
        if(!is.null(gene.list)){
            cat(sprintf("Total %d intersection\n",length(intersect(gene.list$geneID,rownames(ds.res.byTPM)))))
            vv.oneRun <- t(sapply(e.TH,function(tt){
                t(apply(ds.res.byTPM[intersect(gene.list$geneID,rownames(ds.res.byTPM)),c(-1,-2),drop=F], 2,function(x){ sum(x>tt) } ))
            }))
            rownames(vv.oneRun) <- e.TH
        }else{
            vv.oneRun <- t(sapply(e.TH,function(tt){
                t(apply(ds.res.byTPM[,c(-1,-2),drop=F],2,function(x){ sum(x>tt) } ))
            }))
            rownames(vv.oneRun) <- e.TH
        }
        colnames(vv.oneRun) <- colnames(ds.res.byTPM)[c(-1,-2)]
        return(vv.oneRun)
        },.progress = "none",.parallel=T)
    cat(sprintf("[%s] end downsample & count reads:\n", Sys.time()))
    
    save(detected.gn,file=sprintf("%s.RData",out.prefix))
    #save.image(file=sprintf("%s.RData",out.prefix))
}


auto.colSet <- function(n=2,name="Set1")
{
    suppressPackageStartupMessages(require("RColorBrewer"))
    if(n<=9){ 
        ret <- brewer.pal(max(n,3),name)[seq_len(n)] 
    }else{ 
        ret <- colorRampPalette(brewer.pal(9,name))(n) 
    } 
    return(ret)
}

rname <- rownames(detected.gn[[1]])
detected.gn.a <- array(unlist(detected.gn),dim=c(dim(detected.gn[[1]]),length(detected.gn)))
detected.gn.l <- list()
for(i in seq_along(rname)){
    detected.gn.l[[rname[i]]] <- t(detected.gn.a[i,,])
    colnames(detected.gn.l[[rname[i]]]) <- colnames(detected.gn[[1]])
}
detected.gn <- detected.gn.l

#ths <- names(detected.gn)
####e.TH <- c(0,0.1,0.5,1,2,3,5,10,30,50,100)
ths <- c("0","0.5","1","3","10","30","50")
#colSet <- rev(auto.colSet(length(ths)+10,name="Blues"))
colSet <- rev(auto.colSet(length(ths),name="Paired"))
pdf(sprintf("%s.ds.pdf",out.prefix),width=10,height=8)
par(mar=c(10,8,4,10),cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
y <- apply(detected.gn[[1]],2,mean)
yrange <- pretty(1:max(detected.gn[[1]]))
if(args.c){
    x <- as.integer(colnames(detected.gn[[1]]))
    plot(x, y, type="n",xaxt="n",yaxt="s",xlab="",ylab=args.ylab,main=paste0("",sample.id),
         pch=20,lwd=1.5,ylim=c(0,yrange[length(yrange)]))
    axis(1,x,labels = F)
    #staxlab(1,at = setdiff(as.integer(colnames(detected.gn[[1]])),c(150000,200000,250000)),
    #        labels=sprintf("%dK",setdiff(as.integer(colnames(detected.gn[[1]])),c(150000,200000,250000))/1000),
    #        srt=if(length(x) > 5) 90 else 0, cex=1.5,adj=1,top.line=2)
    staxlab(1,at = setdiff(as.integer(colnames(detected.gn[[1]])),c(150000,200000,250000)),
            labels=setdiff(as.integer(colnames(detected.gn[[1]])),c(150000,200000,250000)),
            srt=if(length(x) > 5) 90 else 0, cex=1.5,adj=1,top.line=2)
    i <- 0
    for(nn in ths)
    {
        i <- i+1
        x <- as.integer(colnames(detected.gn[[nn]]))
        y <- apply(detected.gn[[nn]],2,mean)
        y.sd <- apply(detected.gn[[nn]],2,sd)
        lines(x,y,type=ifelse(args.m>1,"l","b"),pch=20,lwd=3,col=colSet[i])
        arrows(x0 = x, y0 = y-y.sd, x1 = x, y1 = y+y.sd,
               length = 0.03,code = 3,angle = 90,lwd = 3,col=colSet[i])
    }
    mtext(sprintf("Read Pair Used",nbanalyzedreads),side=1,line=7.0,cex=2.0)
    ###mtext(sprintf("Reads Used\n(n=%d)",nbanalyzedreads),side=1,line=4.5,cex=2)
}else{
    x <- seq_len(ncol(detected.gn[[1]]))
    plot(x, y, type="n",xaxt = "n",xlab="",ylab=args.ylab,main=paste0("",sample.id),
         pch=20,lwd=1.5,ylim=c(0,yrange[length(yrange)]))
    axis(1,x,colnames(detected.gn[[1]]))
    i <- 5
    for(nn in names(ths))
    {
        i <- i+1
        x <- seq_len(ncol(detected.gn[[nn]]))
        y <- apply(detected.gn[[nn]],2,mean)
        y.sd <- apply(detected.gn[[nn]],2,sd)
        lines(x,y,type=ifelse(args.m>1,"l","b"),pch=20,lwd=1.5,col=colSet[i])
        arrows(x0 = x, y0 = y-y.sd, x1 = x, y1 = y+y.sd,
               length = 0.02,code = 3,angle = 90,lwd = 1.1,col=colSet[i])
    }
    mtext(sprintf("Percentage Of Total Reads\n(n=%d)",nbanalyzedreads),side=1,line=4.5,cex=1.5)
}
legend("right",legend=sprintf("> %s",ths),cex=1.5,
       col=colSet[seq_len(length(ths))+0],lwd=3,pt.cex=3,pch=20,bty="n",inset = -0.25,xpd = T)
#legend("right",legend=sprintf("> %s",names(detected.gn)),cex=2,
#       col=colSet[seq_len(length(detected.gn))+5],lwd=1.5,pch=20,bty="n",inset = -0.25,xpd = T)
dev.off()



