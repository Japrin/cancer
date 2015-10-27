#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<3)
{
	cat("usage: plot_baf_luoxi.R <baf file> <output dir> <sample id>\n")
	q()
}
baf.file <-args[1]
outDir <-args[2]
sampleID <-args[3]
#baf.file <- "S1-1-4.luoxi.baf.txt"
#outDir <- "."
#sampleID <- "S1-1-4"

### read input
LOHMarker <- read.table(baf.file,header=T,check.names=F,stringsAsFactors=F)
names(LOHMarker) <- c("chr","SNP_loc","tumor_BAF")

### chromosome info
chr_list = c(1:22,'X','Y')
chr_length = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
names(chr_length)=chr_list

### translate x coordinate
chr_line_x = c(0)
chr_main_x = c()
x_limit = 0
baf_marker_x_y =c ()
baf_seg_x_y = c()

for(i in 1:length(chr_list))
{
	sChr <- chr_list[i]
	# BAF (marker)
	idx <- which(LOHMarker$chr==sChr)
	t_data<-LOHMarker[idx,]
	baf_marker_x_y=rbind(baf_marker_x_y,data.frame(x=t_data$SNP_loc+x_limit,y=t_data$tumor_BAF,n=t_data$chr,stringsAsFactors=F))
	#
	chr_main_x = c( chr_main_x , mean(c(0,chr_length[sChr]))+x_limit )
	x_limit = x_limit + chr_length[sChr]
	chr_line_x = c( chr_line_x , x_limit )
}

### plot
png(filename=paste(outDir,"/",sampleID,".whole.png",sep=""),width=2500,height=500)
op <- par(mar=c(2,6,6,11)+0.1)
strMain <- paste(sampleID," (#LOH marker=",nrow(LOHMarker),")",sep="")
plot(baf_marker_x_y$x,baf_marker_x_y$y,xlim=c(0,x_limit),ylim = c(0,1),type="p", xlab='',ylab='BAF',xaxt='n', xaxs='i',col="black",pch=16,cex=0.4,bty='n',cex.lab=2.5,cex.axis=2, main=strMain, cex.main=2.5)
abline( v = chr_line_x ,  lty=5,   col = 'gray44')  
abline( h = seq(0,1,0.2), lty = 5 ,col = 'gray44')
text( x = chr_main_x , y = 1.02 , labels = sub("chr","",chr_list), pos = 3 , cex=2, xpd=T)
dev.off()
