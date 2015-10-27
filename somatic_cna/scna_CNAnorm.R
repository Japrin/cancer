#!/usr/local/bin/Rscript

args <- commandArgs(TRUE)
for(i in 1:length(args))
{
	eval(parse(text=args[[i]]))
}

library(CNAnorm)
data(gPar)

##infile
##outfile

setwd(dirname(outfile))

a<-read.table(infile,header=T)

CN <- dataFrame2object(a)

toSkip <- c("chrY", "chrM")
CN <- gcNorm(CN, exclude = toSkip)
CN <- addSmooth(CN, lambda = 7)
CN <- peakPloidy(CN, exclude = toSkip)
png("plotPeaks.png",width=800,height=600)
plotPeaks(CN, special1 = 'chrX', special2 = 'chrY')
dev.off()
CN.default <- validation(CN)
CN <- addDNACopy(CN)
CN <- discreteNorm(CN)


exportTable(CN, file = outfile, show = 'ploidy')
png("all.cnv.png",width=2000,height=800)
plotGenome(CN, superimpose = 'DNACopy', gPar = gPar, colorful = TRUE)
dev.off()

toPlot <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
for(i in 1:length(toPlot))
{
	subSet <- chrs(CN) %in% toPlot[i]
	png(paste(toPlot[i],".cnv.png",sep=""),width=800,height=600)
	plotGenome(CN[subSet], superimpose = 'DNACopy', gPar = gPar, colorful = TRUE)
	dev.off()
}
save.image(file="CNAnorm.RData")
