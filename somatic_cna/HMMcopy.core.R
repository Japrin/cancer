#!/PROJ/GR/share/Software/R/bin/Rscript

args<-commandArgs(T)
if(length(args)<6)
{
	cat("usage: HMMcopy.core.R <normal> <tumor> <gc> <map> <output prefix> <sampleID>\n")
	q()
}
normal_wig<-args[1]
tumor_wig<-args[2]
gc_wig<-args[3]
map_wig<-args[4]
o_prefix<-args[5]
sampleID<-args[6]

library(HMMcopy)

outputBin<-function (correctOutput, segmented_copy, file, sample = "R")
{
	width <- start(correctOutput)[2] - start(correctOutput)[1]
	out <- data.frame(sample = sample, chr = space(correctOutput), start = start(correctOutput) - 1, end = start(correctOutput) + width - 1,gc=correctOutput$gc, map=correctOutput$map, cor.gc=correctOutput$cor.gc, cor.map=correctOutput$cor.map, copy = round(correctOutput$copy, digits = 6),HMMstat = segmented_copy$state )
	write.table(format(out, format = "f", trim = TRUE, drop0trailing = TRUE), file = file, quote = FALSE, row.names = FALSE, sep = "\t")
}

# correction
tum_uncorrected_reads <- wigsToRangedData(tumor_wig, gc_wig, map_wig)
norm_uncorrected_reads <- wigsToRangedData(normal_wig, gc_wig, map_wig)
tum_corrected_copy <- correctReadcount(tum_uncorrected_reads)
norm_corrected_copy <- correctReadcount(norm_uncorrected_reads)
tum_corrected_copy$copy <- tum_corrected_copy$copy - norm_corrected_copy$copy
# Segmenting
param <- HMMsegment(tum_corrected_copy, getparam = TRUE) # retrieve converged parameters via EM
param_default<-param
param$e <- 0.999999999999999
param$strength <- 1e+30
#param$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
#param$m <- param$mu
segmented_copy <- HMMsegment(tum_corrected_copy, param) # perform segmentation via Viterbi
# Export txt
outputBin(tum_corrected_copy, segmented_copy, file = paste(o_prefix,".bin.txt",sep=""),sample=sampleID)
write.table(segmented_copy$segs,file=paste(o_prefix,".segment.txt",sep=""),quote=F,sep="\t",row.names=F)
# Plot
plotBias(tum_corrected_copy)
plotCorrection(tum_corrected_copy,pch=20,col="gray")
plotSegments(tum_corrected_copy, segmented_copy,pch=20)
plotParam(segmented_copy,param)

q(save="yes")
