#!/usr/bin/env Rscript

args = commandArgs(TRUE)

if( length(args)<3 ){
   message("Usage: plot.CNV.WGS.freec.R <freec output> <png/pdf file> <bin size>")
   q()
}

cnv_file<-args[1]
pic_file<-args[2]
bin_size<-as.numeric(args[3])
print(bin_size)
chr_list = c(1:22,'X')
chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560 )
names(chr_length)=chr_list
chr_col = c(rep(c("blue","red"),11),"blue")
y_limit=c(0,4)

# Chromosome      Start   Ratio   MedianRatio     CopyNumber
# 1       4000001 0.917577        0.962999        2
# 1       5000001 0.963934        0.962999        2
cnv_data<-read.table(cnv_file,sep="\t",header=T)

chr_line_x = c(0)
chr_main_x = c()
x_limit = 0

win_logR_x_y =c ()

for(chr in 1:length(chr_list))
{
	t_data<-cnv_data[cnv_data$Chromosome==chr_list[chr],,drop=F]

	# bin
	x<-t_data$Start+x_limit
	y1<-t_data$Ratio
	y2<-t_data$MedianRatio
	y3<-t_data$CopyNumber
	n<-t_data$Chromosome
	y1[y1>=y_limit[2]]<-y_limit[2]
	y2[y2>=y_limit[2]]<-y_limit[2]

	win_logR_x_y=rbind(win_logR_x_y,data.frame(x=x,y1=y1,y2=y2,y3=y3,n=n))
	
	#
	chr_main_x = c( chr_main_x , mean(c(0,chr_length[chr]))+x_limit )
	x_limit = x_limit + chr_length[chr]
	chr_line_x = c( chr_line_x , x_limit )
}
if(regexpr(".png$",pic_file)!=-1)
{
	png(filename=pic_file , width = 2300, height = 400 )
}else
{
	pdf(pic_file,width=23,height=4)
}
par(mar=c(2,6,2,2))
plot( win_logR_x_y[,1],win_logR_x_y[,2],xlim=c(0,x_limit),ylim = c(y_limit[1],y_limit[2]*1.05),type="n", xlab='',ylab='Ratio',xaxt='n', xaxs='i',col="black",pch=16,cex=0.4,bty='n',cex.lab=2.5,cex.axis=2)
#plot( win_logR_x_y[,1],win_logR_x_y[,2],xlim=c(0,x_limit),ylim = c(y_limit[1],y_limit[2]*1.05), xlab='',ylab='Ratio',xaxt='n',col="black",pch=16,cex=0.4,bty='n',cex.lab=2.5,cex.axis=2)
#segments(win_logR_x_y[,1],win_logR_x_y[,3],win_logR_x_y[,1]+bin_size,win_logR_x_y[,3],col="red",lwd=3)
for(chr in 1:length(chr_list))
{	
	tt<-win_logR_x_y[win_logR_x_y$n==chr_list[chr],]
	points(tt[,1],tt[,2],xlim=c(0,x_limit),ylim = c(y_limit[1],y_limit[2]*1.05),col=chr_col[chr],pch=16,cex=0.4,bty='n',cex.lab=2.5,cex.axis=2)
	segments(tt[,1],tt[,3],tt[,1]+bin_size,tt[,3],col="darkgreen",lwd=5)
}

### plot chr_line and add chr_main
chr_num = length( chr_list )
abline( v = chr_line_x ,  lty=5,   col = 'gray44')  
abline( h = seq(y_limit[1],y_limit[2],0.5), lty = 5 ,col = 'gray')
text( x = chr_main_x , y = rep(y_limit[2]*1.04,23) , labels = chr_list , pos = 1 , cex=2)

dev.off()

message("CNV finished")
