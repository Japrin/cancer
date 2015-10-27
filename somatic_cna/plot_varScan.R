#!/usr/local/bin/Rscript

args<-commandArgs(T)

if(length(args)<2)
{
	cat("usage: plot_varScan.R <infile> <outfile> <sampleID>\n")
	q()
}

infile<-args[1]
outfile<-args[2]
if(length(args)>2)
{
	sampleID<-args[3]
}else
{
	sampleID<-"Sample.1"
}

cn <- read.table(infile,header=T)

library(DNAcopy)
CNA.object <-CNA( genomdat = cn$adjusted_log_ratio, chrom = ordered(cn$chrom,levels=c(1:22,"X","Y")), maploc = cn$chr_start, data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=2, min.width=2, undo.splits="sdundo")
#segs2 = segs$output
#write.table(segs2[,2:6], file=outfile, row.names=F, col.names=T, quote=F, sep="\t")
segs_p<-segments.p(segs)
write.table(segs_p, file=outfile, row.names=F, col.names=T, quote=F, sep="\t")

maxLevelToPlot<-3

for (i in c(1:22,'X','Y')) 
{
	png(filename=paste(outfile,".chr",i,".png",sep=""),width=800,height=600,units = "px", pointsize = 20, bg = "white", res = NA)
	tt <- which(cn$chrom==paste("",i,sep=""))
	#print(length(tt))
	if (length(tt)>0) 
	{
		plot(cn$chr_start[tt],cn$adjusted_log_ratio[tt],ylim = c(-maxLevelToPlot,maxLevelToPlot),main=sampleID,xlab = paste ("Position, chr",i),ylab=expression(paste("Log",scriptstyle(2)," Ratio",sep="")),pch = ".",col = "grey") ### colors()[88]

		tt <- which(cn$chrom==paste("",i,sep="") & cn$adjusted_log_ratio >=  maxLevelToPlot)
		points(cn$chr_start[tt],rep(maxLevelToPlot,length(cn$chr_start[tt])),pch = ".",col = "red",cex=4)
		
		tt <- which(cn$chrom==paste("",i,sep="") & cn$adjusted_log_ratio <= -maxLevelToPlot)
		points(cn$chr_start[tt],rep(-maxLevelToPlot,length(cn$chr_start[tt])),pch = ".",col = "blue",cex=4)
		
		#chrom   loc.start       loc.end num.mark        seg.mean
		tt <- which(segs_p$chrom==paste("",i,sep=""))
		ttt <-segs_p[tt,]
		for(i in 1:(dim(ttt)[1]))
		{
			x<-c(ttt$loc.start[i],ttt$loc.end[i])
			y<-c(ttt$seg.mean[i],ttt$seg.mean[i])
			if(ttt$seg.mean[i] > 0.3) 
			{
				color<-"red" 
			}
			else if(ttt$seg.mean[i] < -0.3)
			{
				color<-"blue"
			}
			else
			{
				color<-"black"
			}
			lines(x,y, col = color,lwd=4)
		}

	}
	dev.off()
}
