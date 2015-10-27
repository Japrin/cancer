#.libPaths(new=c("/PUBLIC/software/public/R/v3.0.3/lib64/R/library"))
library(getopt)
library(DNAcopy)
library(ExomeCNV)
suppressMessages(library(doMC))

spec = matrix( c(
   "normal_depth" , "n" , 1 , "character" , "normol_sample gatk_depth_file" ,
   "tumor_depth" , "t" , 1 , "character" , "tumor_sample gatk_depth_file" ,
   "tumor_ID" , "i" , 1 , "character" , "tumor sample id" ,
   "read_length" , "l" , 1 , "integer" , "PE: read_mean_length" ,
   "admixture_rate" , "r" , 1 , "double" , "sample admixture rate" ,
   "output_dir" , "o" , 1 , "character" , "output_dir" ,
    "help"        , "h" , 0 , "logical"   , "help"
   ), ncol=5 , byrow=TRUE )

opt = getopt(spec , opt = commandArgs(TRUE))

if( !is.null(opt$help) ) {
   cat(getopt(spec , usage=TRUE))
   q(status=1)
}

if( is.null(opt$normal_depth) || is.null(opt$tumor_depth) || is.null(opt$tumor_ID) 
   || is.null(opt$read_length) || is.null(opt$admixture_rate) || is.null(opt$output_dir) ){
   cat(getopt(spec , usage=TRUE))
   q()
}

normal_depth = opt$normal_depth
tumor_depth = opt$tumor_depth
tumor_ID = opt$tumor_ID
read_length = opt$read_length
admixture_rate = opt$admixture_rate
output_dir = opt$output_dir

### read depth
normal = read.coverage.gatk(normal_depth)
tumor = read.coverage.gatk(tumor_depth)

n_sum = sum(as.numeric(normal$coverage), na.rm=TRUE)
t_sum = sum(as.numeric(tumor$coverage), na.rm=TRUE)
mean_sum = ( n_sum + t_sum ) / 2 

normal$coverage = ( normal$coverage / n_sum ) * mean_sum
normal$average.coverage = normal$coverage / normal$targeted.base

tumor$coverage = ( tumor$coverage / t_sum ) * mean_sum
tumor$average.coverage = tumor$coverage / tumor$targeted.base

chr.list = unique( normal$chr )

### calculate ratio
logR=calculate.logR(normal,tumor)

### exome CNV
registerDoMC()
eCNV<- foreach(i=1:length(chr.list) , .combine=rbind , .inorder=TRUE) %dopar% {
   idx=which(normal$chr==chr.list[i])
   ecnv=classify.eCNV(normal=normal[idx,], tumor=tumor[idx,], logR=logR[idx], 
      min.spec=0.9999, min.sens=0.9999, option="spec", c=admixture_rate, l=read_length)
}
eCNV_file = paste(sep="",output_dir,'/',tumor_ID,'.eCNV.txt')
write.table(eCNV, file=eCNV_file, quote=FALSE, sep="\t", row.names=FALSE)

### segment CNV
set.seed(50)
cnv=multi.CNV.analyze( normal, tumor, logR=logR, all.cnv.ls=list(eCNV), coverage.cutoff=10, 
   min.spec=0.9995, min.sens=0.9995, option="auc", c=admixture_rate, l=read_length, 
   sdundo=c(1,2), alpha=c(0.05, 0.01) )

# cnv=multi.CNV.analyze( normal, tumor, logR=logR, coverage.cutoff=10,
  #  min.spec=0.999, min.sens=0.999, option="auc", c=admixture_rate, l=read_length, 
  #  sdundo=c(1,2), alpha=c(0.05, 0.01) )

### plot
# temp = paste( sep="", output_dir, '/', tumor_ID, '.eCNV.png' )
# png( temp , width=2000 , height=1200 )
# do.plot.eCNV(eCNV, lim.quantile=0.995, style="idx", line.plot=F)
# dev.off()
# temp = paste( sep="", output_dir, '/', tumor_ID, '.CNV.png' )
# png( temp , width=2000 , height=1200 )
# do.plot.eCNV(cnv, lim.quantile=0.995, style="bp", bg.cnv=e.cnv, line.plot=T)
# dev.off()

write.output(eCNV, cnv, paste(sep="" , output_dir , "/" , tumor_ID), 3)

# `%+%` <- function(x,y) paste(x,y,sep="")
# index = which(!is.nan(cnv$logR) & !cnv$logR %in% c(-Inf,Inf) & cnv$num.mark>=4 & cnv$copy.number!=2)
# write.table(cnv[index ,c("chr","probe_start","probe_end","logR")], file=paste(sep="" , output_dir , "/" , tumor_ID) %+% ".CNV.CBS.lrr.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
