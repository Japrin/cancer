#!/usr/local/bin/Rscript


args<-commandArgs(T)
if(length(args)<4)
{
	cat("usage: plot_ExomeCNV.R <sampleID> <exon.lrr.txt> <cnv.targetN.txt> <outDir> [eLOH.loh.txt] [ExomeCNV.LOH]\n")
	q()
}

sampleID<-args[1]
binFile<-args[2]
segFile<-args[3]
outDir<-args[4]
binTable<-read.table(binFile,header=F)
segTable<-read.table(segFile,header=T,sep="\t")

names(binTable)<-c("chr","beg","end","logR")
binTable$ratio <- 2^(binTable$logR)
#names(segTable)<-c("chr","beg","end","logR")
#cnv.targetN.txt header:
#chr     probe_start     probe_end       coverage        targeted.base   sequenced.base  copy.number     logR    ratio   spec    sens    average.coverage	target.number
#eLOH.loh.txt header:
#chr1    876499  876500  17      26      4       18      FALSE   0.00653470996235409
#ExomeCNV.LOH:
#chr     beg     end     normal_t        normal_v        tumor_t tumor_v isLOH   pValue  target.number   target.LOH.number       mean.mirror.baf mean.LOH.mirror.baf
#chr1    888639  898324  1408    670     1727    559     TRUE    7.83355053541526e-07    8       8       0.787871036257575       0.787871036257575

if(length(args)>=6)
{
	LOHMarkerFile<-args[5]
	LOHSegFile<-args[6]
	LOHMarker<-read.table(LOHMarkerFile,header=F)
	LOHSeg<-read.table(LOHSegFile,header=T)
	names(LOHMarker)<-c("chr","beg","end","normal_v","normal_t","tumor_v","tumor_t","isLOH","pValue")
	#names(LOHSeg)<-c("chr","beg","end","normal_t","normal_v","tumor_t","tumor_v","isLOH","pValue","target.number")
}

chr_list = paste('chr',c(1:22,'X','Y'),sep="")
chr_length = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
names(chr_length)=chr_list


maxLevelToPlot<-3
idx<-which(binTable$ratio>=maxLevelToPlot)
binTable$ratio[idx]<-maxLevelToPlot

## plot by chromoosme
for (i in c(1:22,'X','Y')) 
{
	if(length(args)>=6)
	{
		png(filename=paste(outDir,"/",sampleID,".chr",i,".png",sep=""),width=1400,height=1500,units = "px", pointsize = 20, bg = "white", res = NA)
		op <- par(mfrow = c(3,1),mar=c(5,6,4,11)+0.1)
	}else
	{
		png(filename=paste(outDir,"/",sampleID,".chr",i,".png",sep=""),width=1000,height=600,units = "px", pointsize = 20, bg = "white", res = NA)
		op <- par(mar=c(5,5,4,2)+0.1)
	}
	sChr=paste("chr",i,sep="")
	## cnv plot
	tt <- which(binTable$chr==sChr)
	if(length(tt)>0)
	{
		if(length(args)>=6)
		{
			plot(binTable$beg[tt],binTable$ratio[tt],xlim=c(1,chr_length[sChr]),ylim = c(0,maxLevelToPlot*1.05),main=paste(sampleID," (",sChr,")",sep=""),xlab = "Genomic Position",ylab = "Ratio Profile",pch = 16,col = "grey", cex=0.6, cex.lab=2, cex.main=2.5, cex.axis=1.5)
		}else
		{
			plot(binTable$beg[tt],binTable$ratio[tt],xlim=c(1,chr_length[sChr]),ylim = c(0,maxLevelToPlot*1.05),main=paste(sampleID," (",sChr,")",sep=""),xlab = "Genomic Position",ylab = "Ratio Profile",pch = 16,col = "grey", cex=0.6, cex.lab=1.25, cex.main=1.5, cex.axis=1.0)
		}
		abline( h = seq(0,maxLevelToPlot,0.5), lty = 5 ,col = 'gray44')
		
		tt <- which(binTable$chr==paste("chr",i,sep="") & binTable$ratio>=maxLevelToPlot)
		points(binTable$beg[tt],rep(maxLevelToPlot,length(binTable$beg[tt])),pch = 6,col = "red",cex=0.5)
		 
		#tt <- which(binTable$chr==paste("chr",i,sep="") & binTable$logR<=-maxLevelToPlot)
		#points(binTable$beg[tt],rep(-maxLevelToPlot,length(binTable$beg[tt])),pch = 2,col = "blue",cex=0.5)
		
		tt <- which(segTable$chr==paste("chr",i,sep=""))
		ttt <-segTable[tt,]
		for(j in 1:(dim(ttt)[1]))
		{
			x<-c(ttt$probe_start[j],ttt$probe_end[j])
			y<-c(ttt$ratio[j],ttt$ratio[j])
			color <- "grey"
			if(is.na(ttt$copy.number[j]) || ttt$copy.number[j] == 2)
			{
				color<-"black"
			}else if(ttt$copy.number[j] > 2 && ttt$target.number[j] >= 10 ) 
			{
				color<-"red" 
			}
			else if(ttt$copy.number[j] < 2 && ttt$target.number[j] >= 10 )
			{
				color<-"blue"
			}
			lines(x,y, col = color,lwd=4)
		}
		legend(x="topright",legend=paste("CN=",c(1:3),sep=""),col=c("blue","green","red"), pch=20, xpd=NA, inset=c(-0.12,0.00),cex=2)
	}
	## LOH plot
	if(length(args)>=6)
	{
		tt <- which(LOHMarker$chr==sChr)
		ss <- which(LOHSeg$chr==sChr)
		ssLOH <- which(LOHSeg$chr==sChr & LOHSeg$isLOH & LOHSeg$target.number >=5)
		#ssLOH <- which(LOHSeg$chr==sChr & LOHSeg$isLOH )
		ssNoLOH <- which(LOHSeg$chr==sChr & (!LOHSeg$isLOH ) )
		#ssNoLOH <- which(LOHSeg$chr==sChr & (!LOHSeg$isLOH | LOHSeg$target.number <5) )
		#print(LOHSeg[ssLOH,])
		plot(LOHMarker$beg[tt],LOHMarker$tumor_v[tt]/LOHMarker$tumor_t[tt],xlim=c(1,chr_length[sChr]),ylim = c(0,1),main="",xlab = "Genomic Position",ylab = "BAF Profile (Tumor)",pch = 16,col = "grey", cex=0.6, cex.lab=2, cex.main=2.5, cex.axis=1.5)
		legend(x="topright",legend=c("LOH ","HET "),col=c("red","green"), pch=20, xpd=NA, inset=c(-0.12,0.00),cex=2)
		abline( h = seq(0,1,0.2), lty = 5 ,col = 'gray44')
		if(length(ssLOH)>0)
		{
			segY <- LOHSeg$mean.mirror.baf[ssLOH]
			segments(LOHSeg$beg[ssLOH],segY,LOHSeg$end[ssLOH],segY,col="orange",lwd=5)
			segments(LOHSeg$beg[ssLOH],1-segY,LOHSeg$end[ssLOH],1-segY,col="orange",lwd=5)
		}
		if(length(ssNoLOH)>0)
		{
			segY <- LOHSeg$mean.mirror.baf[ssNoLOH]
			#segments(LOHSeg$beg[ssNoLOH],segY,LOHSeg$end[ssNoLOH],segY,col="black",lwd=5)
			#segments(LOHSeg$beg[ssNoLOH],1-segY,LOHSeg$end[ssNoLOH],1-segY,col="black",lwd=5)
		}
		plot(LOHMarker$beg[tt],LOHMarker$normal_v[tt]/LOHMarker$normal_t[tt],xlim=c(1,chr_length[sChr]),ylim = c(0,1),main="",xlab = "Genomic Position",ylab = "BAF Profile (Normal)",pch = 16,col = "grey", cex=0.6, cex.lab=2, cex.main=2.5, cex.axis=1.5)
		abline( h = seq(0,1,0.2), lty = 5 ,col = 'gray44')
		segY <- LOHSeg$mean.mirror.baf.normal[ss]
		#segments(LOHSeg$beg[ss],segY,LOHSeg$end[ss],segY,col="black",lwd=5)
		#segments(LOHSeg$beg[ss],1-segY,LOHSeg$end[ss],1-segY,col="black",lwd=5)
	}
	dev.off()
}

## plot whole genome
chr_line_x = c(0)
chr_main_x = c()
x_limit = 0

win_x_y =c ()
seg_x_y =c ()
baf_marker_x_y =c ()
baf_seg_x_y = c()

for(i in 1:length(chr_list))
{
	sChr <- chr_list[i]
	# bin
	idx <- which(binTable$chr==sChr)
	t_data<-binTable[idx,]
	win_x_y=rbind(win_x_y,data.frame(x=t_data$beg+x_limit,y1=t_data$ratio,n=t_data$chr))
	# seg
	idx <- which(segTable$chr==sChr)
	t_data<-segTable[idx,]
	seg_x_y=rbind(seg_x_y,data.frame(x1=t_data$probe_start+x_limit,x2=t_data$probe_end+x_limit,y1=t_data$ratio,y2=t_data$ratio,n=t_data$chr,copy_number=t_data$copy.number,target_number=t_data$target.number))
	if(length(args)>=6)
	{
	    # BAF (marker)
	    idx <- which(LOHMarker$chr==sChr)
	    t_data<-LOHMarker[idx,]
	    baf_marker_x_y=rbind(baf_marker_x_y,data.frame(x=t_data$beg+x_limit,y=t_data$tumor_v/t_data$tumor_t,n=t_data$chr))
	    # BAF (seg)
	    idx <- which(LOHSeg$chr==sChr)
	    t_data<-LOHSeg[idx,]
	    baf_seg_x_y=rbind(baf_seg_x_y,data.frame(x1=t_data$beg+x_limit,x2=t_data$end+x_limit,y1=t_data$mean.mirror.baf,y2=1-t_data$mean.mirror.baf,isLOH=t_data$isLOH,target.number=t_data$target.number,n=t_data$chr))
	}
	#
	chr_main_x = c( chr_main_x , mean(c(0,chr_length[sChr]))+x_limit )
	x_limit = x_limit + chr_length[sChr]
	chr_line_x = c( chr_line_x , x_limit )
}

if(length(args)>=6)
{
	png(filename=paste(outDir,"/",sampleID,".whole.png",sep=""),width=2500,height=900)
	layout(matrix(c(1,2)),heights=c(0.6,0.4))
	op <- par(mar=c(2,6,6,11)+0.1)
	strMain <- paste(sampleID," (#CNV marker=",nrow(binTable),", #LOH marker=",nrow(LOHMarker),")",sep="")
}else
{
	png(filename=paste(outDir,"/",sampleID,".whole.png",sep=""),width=2500,height=500)
	op <- par(mar=c(2,6,6,11)+0.1)
	strMain <- paste(sampleID," (#CNV marker=",nrow(binTable),")",sep="")
}

plot( win_x_y$x,win_x_y$y1,xlim=c(0,x_limit),ylim = c(0,maxLevelToPlot*1.05),type="n", xlab='',ylab='Ratio',xaxt='n', xaxs='i',col="black",pch=16,cex=0.4,bty='n',cex.lab=2.5,cex.axis=2, main=strMain, cex.main=2.5)
legend(x="topright",legend=paste("CN=",c(1:3),sep=""),col=c("blue","green","red"), pch=20, xpd=NA, inset=c(-0.055,0.00),cex=2)

for(i in 1:length(chr_list))
{
	sChr=chr_list[i]
	tt<-win_x_y[win_x_y$n==sChr,]
	points(tt$x,tt$y1,col="gray",pch=16,cex=0.4)
	tt<-seg_x_y[seg_x_y$n==sChr & seg_x_y$copy_number==2,]
	segments(tt$x1,tt$y1,tt$x2,tt$y2,col="darkgreen",lwd=5)
	tt<-seg_x_y[seg_x_y$n==sChr & seg_x_y$copy_number>2 & seg_x_y$target_number>=10,]
	segments(tt$x1,tt$y1,tt$x2,tt$y2,col="red",lwd=5)
	tt<-seg_x_y[seg_x_y$n==sChr & seg_x_y$copy_number<2 & seg_x_y$target_number>=10,]
	segments(tt$x1,tt$y1,tt$x2,tt$y2,col="blue",lwd=5)
}

### plot chr_line and add chr_main
chr_num = length( chr_list )
abline( v = chr_line_x ,  lty=5,   col = 'gray44')  
abline( h = seq(0,maxLevelToPlot,0.5), lty = 5 ,col = 'gray44')
text( x = chr_main_x , y = maxLevelToPlot*1.05 , labels = sub("chr","",chr_list), pos = 3 , cex=2, xpd=T)

if(length(args)>=6)
{
	par(mar=c(2,6,2,11)+0.1)
	plot(baf_marker_x_y$x,baf_marker_x_y$y,xlim=c(0,x_limit),ylim = c(0,1),type="n", xlab='',ylab='BAF',xaxt='n', xaxs='i',col="black",pch=16,cex=0.4,bty='n',cex.lab=2.5,cex.axis=2)
	legend(x="topright",legend=c("LOH ","HET "),col=c("red","green"), pch=20, xpd=NA, inset=c(-0.055,0.00),cex=2)
	for(i in 1:length(chr_list))
	{
	    sChr=chr_list[i]
	    tt<-baf_marker_x_y[baf_marker_x_y$n==sChr,]
	    points(tt$x,tt$y,col="gray",pch=16,cex=0.4)
	    tt<-baf_seg_x_y[baf_seg_x_y$n==sChr & baf_seg_x_y$isLOH & baf_seg_x_y$target.number>=5,]
	    segments(tt$x1,tt$y1,tt$x2,tt$y1,col="orange",lwd=5)
	    segments(tt$x1,tt$y2,tt$x2,tt$y2,col="orange",lwd=5)
	}
	abline( v = chr_line_x ,  lty=5,   col = 'gray44')  
	abline( h = seq(0,1,0.2), lty = 5 ,col = 'gray44')
	#text( x = chr_main_x , y = maxLevelToPlot*1.05 , labels = sub("chr","",chr_list), pos = 3 , cex=2, xpd=T)

}

dev.off()
