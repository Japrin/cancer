#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--in",required=T)
parser$add_argument("-o", "--out",required=T)
parser$add_argument("-a", "--a",default=0.05,type="double",help="downsample range's start point [default %(default)s]")
parser$add_argument("-b", "--b",default=0.05,type="double",help="downsample range's step [default %(default)s]")
parser$add_argument("-p", "--pair", action="store_true", default=TRUE, help="if specified, PE reads other wize SE reads [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose output [default %(default)s]")
parser$add_argument("-s", "--sample", default="SAMPLE", help="sample id [default %(default)s]")
parser$add_argument("-g", "--gfeature", default="/home/zhenglt/work/TCR_chunhong/technological/bin/hg19.knownGene.RData", help="RData contains gene feature [default %(default)s]")
args <- parser$parse_args()

in.bam.file <- args[["in"]]
out.prefix <- args$out
paired_ends <- args$pair
verbose <- args$verbose
sample.id <- args$sample
gfeature.file <- args$gfeature
ds.start <- args$a
ds.step <- args$b

#### For test only
#in.bam.file <- "../phase04/OUT/TTH8-0205/bams/TTH8-0205.analyzed.bam"
#out.prefix <- "./test.geneNumber_saturation"
#paired_ends <- TRUE
#verbose <- TRUE
#sample.id <- "SAMPLE"
#gfeature.file <- "./bin/hg19.knownGene.RData"
##gfeature.file <- "/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/hg19.knownGene.RData"
#ds.start <- 0.2
#ds.step <- 0.2

#### functin definition
suppressPackageStartupMessages(library("HTSeqGenie"))
suppressPackageStartupMessages(library("RColorBrewer"))

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

#### process input data

reads <- readRNASeqEnds(in.bam.file)
#head(reads)
genomic_features <- get(load(gfeature.file))
## consolidate (i.e. pair) read ends
if (paired_ends) reads <- consolidateByRead(reads)
## remove read names to save ~50% of memory space
reads <- unname(reads)
gc()

nbanalyzedreads <- length(reads)
ds <- sample(1:nbanalyzedreads)
#ds.pct <- seq(0.05,1,0.05)
ds.pct <- seq(ds.start,1,ds.step)
ds.n <- as.integer(nbanalyzedreads*ds.pct)

ds.res.byRPKM <- data.frame()
ds.res.byTPM <- data.frame()
ds.res.byCount <- data.frame() 
hh <- sapply(seq_along(ds.n),function(i,n,pct){
    res.ds <- countGenomicFeaturesSimple(features = genomic_features$gene_exonic,reads = reads[ds[1:n[i]]],nbanalyzedreads = n[i])
    if(ncol(ds.res.byRPKM)==0)
    {
	    ds.res.byRPKM <<- res.ds[,c("name","width")]
	    ds.res.byTPM <<- res.ds[,c("name","width")]
	    ds.res.byCount <<- res.ds[,c("name","width")]
    }
    ds.res.byRPKM[,paste0("",pct[i],"")] <<- res.ds$rpkm
    ds.res.byTPM[,paste0("",pct[i],"")] <<- res.ds$tpm
    ds.res.byCount[,paste0("",pct[i],"")] <<- res.ds$count
    
},ds.n,ds.pct)

if(verbose)
{
	write.table(ds.res.byRPKM,paste0(out.prefix,".ds.RPKM.txt"),sep="\t",quote = F,row.names = F)
	write.table(ds.res.byTPM,paste0(out.prefix,".ds.TPM.txt"),sep="\t",quote = F,row.names = F)
	write.table(ds.res.byCount,paste0(out.prefix,".ds.Count.txt"),sep="\t",quote = F,row.names = F)
}

plot.geneNumber.saturation <- function(ds.res,note="RPKM")
{

    out.gnumber <- data.frame()
    for(i in c(0.1,1,2,3,5,10,50,100,500))
    {
        out.gnumber <- rbind(out.gnumber,
                                    data.frame(sample=sample.id,
                                               filter=paste0(note,">=",i),
                                               t(apply(ds.res[,c(-1,-2)],2,function(x,tt){ sum(x>=tt) },i )),
                                               check.names = F))
    }
    write.table(out.gnumber,paste0(out.prefix,".gn.",note,".txt"),sep="\t",quote = F,row.names = F)

    colSet <- brewer.pal(9,"Set1")
    pdf(paste0(out.prefix,".ds.",note,".pdf"),width=10,height=8)
    par(mar=c(8,6,4,10),cex.lab=1.5)
    plot(colnames(out.gnumber[,c(-1,-2)]),
         out.gnumber[1,c(-1,-2)],
         type="n",xlab="",ylab="Number Of Expressed Genes",main=paste0("",sample.id),pch=20,lwd=1.5,
         ylim = c(0,max(out.gnumber[,c(-1,-2)]))
         )
    for(i in 1:nrow(out.gnumber))
    {
        x=out.gnumber[i,]
        lines(colnames(x[c(-1,-2)]),x[c(-1,-2)],type="b",pch=20,lwd=1.5,col=colSet[as.numeric(x[2])])
    }
    legend("right",legend=out.gnumber$filter,col=colSet[as.numeric(out.gnumber$filter)],lwd=1.5,pch=20,bty="n",inset = -0.25,xpd = T)
    mtext(paste0("Percentage Of Total Reads\n(n=",nbanalyzedreads,")"),side=1,line=4.5,cex=1.5)
    dev.off()
}
plot.geneNumber.saturation(ds.res.byRPKM,note="RPKM")
plot.geneNumber.saturation(ds.res.byTPM,note="TPM")
