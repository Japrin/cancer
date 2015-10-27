#!/usr/local/bin/Rscript


args<-commandArgs(T)
if(length(args)<3)
{
	cat("usage: plot_CONTRA.R <bin.txt> <bin.txt.LargeDeletion.txt> <outDir>\n")
	q()
}

binFile<-args[1]
segFile<-args[2]
outDir<-args[3]

binTable<-read.table(binFile,header=T)
segTable<-read.table(segFile,header=T)

maxLevelToPlot<-3
for (i in c(1:length(binTable$Adjusted.Mean.of.LogRatio))) {
	if (!is.na(binTable$Adjusted.Mean.of.LogRatio[i]) && binTable$Adjusted.Mean.of.LogRatio[i]>maxLevelToPlot) {
		binTable$Adjusted.Mean.of.LogRatio[i]=maxLevelToPlot;
	}
}

for (i in c(1:22,'X','Y')) 
{
	#print(paste("Format.chr",i,".png",sep=""))
	png(filename=paste(outDir,"/chr",i,".png",sep=""),width=800,height=600,units = "px", pointsize = 20, bg = "white", res = NA)
	tt <- which(binTable$Chr==paste("chr",i,sep=""))
	#print(length(tt))
	if (length(tt)>0) 
	{
		plot(binTable$OriStCoordinate[tt],binTable$Adjusted.Mean.of.LogRatio[tt],ylim = c(-maxLevelToPlot,maxLevelToPlot),xlab = paste ("position, chr",i),ylab = "log ratio profile",pch = ".",col = "grey") ### colors()[88]
		
		tt <- which(binTable$Chr==paste("chr",i,sep="") & binTable$gain.loss =="gain" & binTable$Adjusted.P.Value<0.1)
		#tt <- which(binTable$Chr==paste("chr",i,sep="") & binTable$gain.loss =="gain" & binTable$P.Value<0.05)
		points(binTable$OriStCoordinate[tt],binTable$Adjusted.Mean.of.LogRatio[tt],pch = ".",col = colors()[136],cex=4)
		
		tt <- which(binTable$Chr==paste("chr",i,sep="") & binTable$Adjusted.Mean.of.LogRatio>=maxLevelToPlot)
		points(binTable$OriStCoordinate[tt],rep(maxLevelToPlot,length(binTable$OriStCoordinate[tt])),pch = ".",col = colors()[136],cex=4)
		 
		tt <- which(binTable$Chr==paste("chr",i,sep="") & binTable$gain.loss =="loss" & binTable$Adjusted.P.Value<0.1)
		#tt <- which(binTable$Chr==paste("chr",i,sep="") & binTable$gain.loss =="loss" & binTable$P.Value<0.05)
		points(binTable$OriStCoordinate[tt],binTable$Adjusted.Mean.of.LogRatio[tt],pch = ".",col = colors()[461],cex=4)
		
		tt <- which(binTable$Chr==paste("chr",i,sep="") & binTable$Adjusted.Mean.of.LogRatio<=-maxLevelToPlot)
		points(binTable$OriStCoordinate[tt],rep(-maxLevelToPlot,length(binTable$OriStCoordinate[tt])),pch = ".",col = colors()[24],cex=4)
		
		
		tt <- which(segTable$Chr==paste("chr",i,sep=""))
		ttt <-segTable[tt,]
		for(i in 1:(dim(ttt)[1]))
		{
			x<-c(ttt$OriStCoordinate[i],ttt$OriEndCoordinate[i])
			y<-c(ttt$CBS.Mean[i],ttt$CBS.Mean[i])
			lines(x,y, col = colors()[463],lwd=4)
		}
	}
	dev.off()
}
