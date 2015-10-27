#!/usr/local/bin/Rscript

args <- commandArgs(T)

if(length(args)<1)
{
	cat("plot_CNA_ratio.R <infile>\n")
	q()
}
dataTable <-read.table(args[1], header=TRUE);

ratio<-data.frame(dataTable)
fa<-1
if(length(args)>1)
{
	fa <- type.convert(args[2])
}
#head(ratio)
maxLevelToPlot <- 3
for (i in c(1:length(ratio$v_window))) {
	if (!is.na(ratio$v_window[i]) && ratio$v_window[i]>maxLevelToPlot) {
		ratio$v_window[i]=maxLevelToPlot;
	}
}

#chr     beg     end     length  v_type  v_winow v_segment
for (i in c(1:22,'X','Y')) 
{
	#print(paste("Format.chr",i,".png",sep=""))
	png(filename=paste("Format.chr",i,".png",sep=""),width=800,height=600,units = "px", pointsize = 20, bg = "white", res = NA)
	tt <- which(ratio$chr==i)
	if (length(tt)>0) 
	{
		plot(ratio$beg[tt],ratio$v_window[tt]*fa,ylim = c(0,maxLevelToPlot*fa),xlab = paste ("position, chr",i),ylab = "normalized copy number profile",pch = ".",col = colors()[88])
		
		tt <- which(ratio$chr==i  & ratio$v_type=="gain" )
		points(ratio$beg[tt],ratio$v_window[tt]*fa,pch = ".",col = colors()[136])
		
		tt <- which(ratio$chr==i  & ratio$v_window>=maxLevelToPlot)	
		#points(ratio$beg[tt],maxLevelToPlot*fa,pch = ".",col = colors()[136],cex=4)
		 
		tt <- which(ratio$chr==i  & ratio$v_type=="loss")
		points(ratio$beg[tt],ratio$v_window[tt]*fa,pch = ".",col = colors()[461])
		
		tt <- which(ratio$chr==i  & ratio$v_window<=-maxLevelToPlot)
		#points(ratio$beg[tt],ratio$v_window[tt]*fa, pch = ".", col = colors()[24],cex=4)
		
		tt <- which(ratio$chr==i)
		#points(ratio$beg[tt],ratio$v_segment[tt]*fa, pch = ".", col = colors()[463],cex=4)
	}
	dev.off()
}
