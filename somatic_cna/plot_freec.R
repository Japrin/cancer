#!/usr/bin/env Rscript


args<-commandArgs(T)
if(length(args)<4)
{
	cat("usage: plot_freec.R <sampleID> <cnv.result> <cnv.result.segment> <outDir> [zygosity.res] [zygosity.res.segment]\n")
	q()
}

sampleID<-args[1]
binFile<-args[2]
segFile<-args[3]
outDir<-args[4]
binTable<-read.table(binFile,header=T,sep="\t")
segTable<-read.table(segFile,header=T,sep="\t")
names(binTable)<-c("chr","beg","end","copy.number","tumor_DOC","control_DOC","ratio","seg_beg","seg_end","seg_ratio")
binTable$chr <- sub("^chr","",binTable$chr,perl = T)
segTable$chr <- sub("^chr","",segTable$chr,perl = T)
col_mat<-col2rgb(binTable$copy.number+1)
binTable$col=rgb(col_mat[1,],col_mat[2,],col_mat[3,],250,maxColorValue=255)

#cnv.result:
#chr     exon_start      exon_end        cnv     tumor_DOC       control_DOC     rationormalized_after_smoothing CNV_start       CNV_end       seg_mean
#1       762097  762270  3       910     735     1.13376066560237        762097  8716449 1.04741451911632
#1       861281  861490  3       129     100     1.05699725331246        762097  8716449 1.04741451911632
#cnv.result.segment:
#chr     beg     end     ratio   copy.number     target.number
#1       762097  8716449 1.04741451911632        3       962
#1       8921238 9006013 1.47166268731018        5       12
#zygosity.res:
#chrom   SNP_loc control_BAF     tumor_BAF       mirrored_BAF    control_doc     tumor_doc       control_vdoc    tumor_vdoc      cn   z	zygosity
#1       876499  0.653846        0.222222        0.222222        26      18      17      4       3       1       LOH
#1       877715  0.45    0.366667        0.366667        20      30      9       11      3       1       LOH
#zygosity.res.segment:
#chr     beg     end     zygosity        mean.mirror.baf target.number
#1       876499  17046613        LOH     0.2345741       330
#1       17248500        17413121        ASCNA   0.454710111111111       9

if(length(args)>=6)
{
	LOHMarkerFile<-args[5]
	LOHSegFile<-args[6]
	LOHMarker<-read.table(LOHMarkerFile,header=T)
	LOHSeg<-read.table(LOHSegFile,header=T)
	LOHMarker$chr <- sub("^chr","",LOHMarker$chr,perl = T)
	LOHSeg$chr <- sub("^chr","",LOHSeg$chr,perl = T)
	col_mat<-col2rgb(LOHMarker$z+1)
	LOHMarker$col=rgb(col_mat[1,],col_mat[2,],col_mat[3,],250,maxColorValue=255)
}

chr_list = c(1:22,'X','Y')
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
		png(filename=paste(outDir,"/",sampleID,".chr",i,".png",sep=""),width=1400,height=600,units = "px", pointsize = 20, bg = "white", res = NA)
		op <- par(mar=c(5,5,4,11)+0.1)
	}
	#sChr=paste("chr",i,sep="")
	sChr=i
	## cnv plot
	tt <- which(binTable$chr==sChr)
	if(length(tt)>0)
	{
		if(length(args)>=6)
		{
			plot(binTable$beg[tt],binTable$ratio[tt],xlim=c(1,chr_length[sChr]),ylim = c(0,maxLevelToPlot*1.05),main=paste(sampleID," (chr",sChr,")",sep=""),xlab = "Genomic Position",ylab = "Ratio Profile",pch = 16,col = binTable$col[tt], cex=0.6, cex.lab=2, cex.main=2.5, cex.axis=1.5)
		}else
		{
			plot(binTable$beg[tt],binTable$ratio[tt],xlim=c(1,chr_length[sChr]),ylim = c(0,maxLevelToPlot*1.05),main=paste(sampleID," (chr",sChr,")",sep=""),xlab = "Genomic Position",ylab = "Ratio Profile",pch = 16,col = binTable$col[tt], cex=0.6, cex.lab=1.25, cex.main=1.5, cex.axis=1.0)
		}
		abline( h = seq(0,maxLevelToPlot,0.5), lty = 5 ,col = 'gray44')
		
		tt <- which(binTable$chr==sChr & binTable$ratio>=maxLevelToPlot)
		points(binTable$beg[tt],rep(maxLevelToPlot,length(binTable$beg[tt])),pch = 6,col = "red",cex=0.5)
		
		tt <- which(segTable$chr==sChr)
		ttt <-segTable[tt,]
	    	segments(ttt$beg,ttt$ratio,ttt$end,ttt$ratio,col="black",lwd=4)
    		
		legend(x="topright",legend=paste("CN=",as.numeric(names(table(binTable$copy.number))),sep=""),col=as.integer(names(table(binTable$copy.number)))+1, pch=20, xpd=NA, inset=c(-0.12,0.00),cex=2)
	}
	## LOH plot
	if(length(args)>=6)
	{
		tt <- which(LOHMarker$chr==sChr)
		ss <- which(LOHSeg$chr==sChr)
		plot(LOHMarker$SNP_loc[tt],LOHMarker$tumor_BAF[tt],xlim=c(1,chr_length[sChr]),ylim = c(0,1),main="",xlab = "Genomic Position",ylab = "BAF Profile (Tumor)",pch = 16,col = LOHMarker$col[tt], cex=0.6, cex.lab=2, cex.main=2.5, cex.axis=1.5)
		abline( h = seq(0,1,0.2), lty = 5 ,col = 'gray44')
		segments(LOHSeg$beg[ss],LOHSeg$mean.mirror.baf[ss],LOHSeg$end[ss],LOHSeg$mean.mirror.baf[ss],col="black",lwd=5)
		segments(LOHSeg$beg[ss],1-LOHSeg$mean.mirror.baf[ss],LOHSeg$end[ss],1-LOHSeg$mean.mirror.baf[ss],col="black",lwd=5)
    		
		legend(x="topright",legend=c("LOH","HET","ASCNA"),col=c(2:4), pch=20, xpd=NA, inset=c(-0.14,0.00),cex=2)
		
		plot(LOHMarker$SNP_loc[tt],LOHMarker$control_BAF[tt],xlim=c(1,chr_length[sChr]),ylim = c(0,1),main="",xlab = "Genomic Position",ylab = "BAF Profile (Normal)",pch = 16,col = 3, cex=0.6, cex.lab=2, cex.main=2.5, cex.axis=1.5)
		abline( h = seq(0,1,0.2), lty = 5 ,col = 'gray44')
	}
	dev.off()
}

# plot whole genome
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
	if(sChr=="Y" && length(idx)<1)
	{
		chr_list = c(1:22,'X')
		break
	}
	
	t_data<-binTable[idx,]
	win_x_y=rbind(win_x_y,data.frame(x=t_data$beg+x_limit,y1=t_data$ratio,n=t_data$chr,col=t_data$col,stringsAsFactors=F))
	# seg
	idx <- which(segTable$chr==sChr)
	#if(length(idx)<1) next
	t_data<-segTable[idx,]
	seg_x_y=rbind(seg_x_y,data.frame(x1=t_data$beg+x_limit,x2=t_data$end+x_limit,y1=t_data$ratio,y2=t_data$ratio,n=t_data$chr,copy_number=t_data$copy_number,target_number=t_data$target_number,stringsAsFactors=F))
	if(length(args)>=6)
	{
	    # BAF (marker)
	    idx <- which(LOHMarker$chr==sChr)
	    #if(length(idx)<1) next
	    t_data<-LOHMarker[idx,]
	    baf_marker_x_y=rbind(baf_marker_x_y,data.frame(x=t_data$SNP_loc+x_limit,y=t_data$tumor_BAF,n=t_data$chr,col=t_data$col,stringsAsFactors=F))
	    # BAF (seg)
	    idx <- which(LOHSeg$chr==sChr)
	    t_data<-LOHSeg[idx,]
	    baf_seg_x_y=rbind(baf_seg_x_y,data.frame(x1=t_data$beg+x_limit,x2=t_data$end+x_limit,y1=t_data$mean.mirror.baf,y2=1-t_data$mean.mirror.baf,target.number=t_data$target.number,n=t_data$chr))
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

## point size
if(nrow(binTable)>1e6) binCex=0.2  else binCex=0.4 
if(nrow(LOHMarker)>1e6) lohCex=0.2 else lohCex=0.4

plot(win_x_y$x,win_x_y$y1,xlim=c(0,x_limit),ylim = c(0,maxLevelToPlot*1.05),type="p", xlab='',ylab='Ratio',xaxt='n', xaxs='i',col=win_x_y$col,pch=16,cex=binCex,bty='n',cex.lab=2.5,cex.axis=2, main=strMain, cex.main=2.5)
#segments(seg_x_y$x1,seg_x_y$y1,seg_x_y$x2,seg_x_y$y2,col="black",lwd=5)

legend(x="topright",legend=paste("CN=",as.numeric(names(table(binTable$copy.number))),sep=""),col=as.integer(names(table(binTable$copy.number)))+1, pch=20, xpd=NA, inset=c(-0.055,0.00),cex=2)
### plot chr_line and add chr_main
chr_num = length( chr_list )
abline( v = chr_line_x ,  lty=5,   col = 'gray44')  
abline( h = seq(0,maxLevelToPlot,0.5), lty = 5 ,col = 'gray44')
text( x = chr_main_x , y = maxLevelToPlot*1.05 , labels = sub("chr","",chr_list), pos = 3 , cex=2, xpd=T)

if(length(args)>=6)
{
	par(mar=c(2,6,2,11)+0.1)
	plot(baf_marker_x_y$x,baf_marker_x_y$y,xlim=c(0,x_limit),ylim = c(0,1),type="p", xlab='',ylab='BAF',xaxt='n', xaxs='i',col=baf_marker_x_y$col,pch=16,cex=lohCex,bty='n',cex.lab=2.5,cex.axis=2)
	#segments(baf_seg_x_y$x1,baf_seg_x_y$y1,baf_seg_x_y$x2,baf_seg_x_y$y1,col="black",lwd=5)
	#segments(baf_seg_x_y$x1,baf_seg_x_y$y2,baf_seg_x_y$x2,baf_seg_x_y$y2,col="black",lwd=5)
	
	abline( v = chr_line_x ,  lty=5,   col = 'gray44')  
	abline( h = seq(0,1,0.2), lty = 5 ,col = 'gray44')

	legend(x="topright",legend=c("LOH","HET","ASCNA"),col=c(2:4), pch=20, xpd=NA, inset=c(-0.065,0.00),cex=2)
}

dev.off()
