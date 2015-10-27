.libPaths(new=c("/WORK/GR/zhanglanying/01software/Rlib"))
suppressMessages(library(doMC))
### For ExomeCNV result
### plot all_length for each chromosome

args = commandArgs(TRUE)
if(length(args) !=3 ){
   message("Usage: /usr/lib64/R/bin/Rscript /WORK/GR/zhanglanying/program/analysis-pipeline/ExomeCNV/exomeCNV.plot.LOH.23.R <.eLOH.txt> <.LOH.txt> <png_file>")
   q()
}
eLOH_txt = args[1]
loh_txt = args[2]
png_file = args[3]
if(length(args)==4){
   num_mark = args[4]
}else{
   num_mark = 4
}


# eLOH_txt = "/WORK/GR/zhanglanying/work/WGS.breastCaner.3/LOH/189NT.eLOH.txt"
# loh_txt = "/WORK/GR/zhanglanying/work/WGS.breastCaner.3/LOH/189NT.LOH.txt"

### chr base_information
chr_list = paste('chr',1:22,sep="")
chr_list = c(chr_list , 'chrX')
chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 
133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560 )
names(chr_length)=chr_list

chr_length_cum = cumsum(chr_length)
chr_length_cum = c(0, chr_length_cum)

#########  Read Data  #########
#########  Read Data  #########
eLOH =  read.table(eLOH_txt, sep="\t", header=TRUE, as.is=TRUE)
eLOH = eLOH[(eLOH$chr != "chrY"),]

if(file.info(loh_txt)$size == 0){
   png( filename=png_file , width = 400, height = 400 )
   plot(1:10,1:10,type="n",xaxt='n',yaxt='n',xlab="",ylab="",bty='n')
   dev.off()
   q()
}


loh = read.table(loh_txt, sep="\t", header=TRUE, as.is=TRUE)
loh = loh[(loh$LOH=="TRUE" & loh$num.mark>=num_mark & loh$chr != "chrY"),]

normal_baf_y = eLOH$normal.baf / eLOH$normal.coverage
tumor_baf_y = cbind(eLOH$tumor.baf / eLOH$tumor.coverage, eLOH$chr)

### prepare plot data ###
### prepare plot data ###

registerDoMC()
### BAF_X
message("BAF_X")
baf_x<- foreach(i=1:length(chr_list) , .combine=rbind , .inorder=TRUE) %dopar% {
   chr_baf_x = eLOH[(eLOH$chr==chr_list[i]),"position"]
   message(paste(sep=" : ", chr_list[i], length(chr_baf_x) ))
   baf_x_chr = c()
   if(length(chr_baf_x) != 0){
      baf_x_chr = cbind( chr_baf_x + chr_length_cum[i], rep(chr_list[i], length(chr_baf_x)) )
   }
   baf_x_chr = baf_x_chr
}

### LOH_X_Y
message("LOH_X_Y")
tumor_loh_x_y <- foreach(n=1:length(chr_list) , .combine=rbind , .inorder=TRUE) %dopar% {
   chr_loh_x = loh[loh$chr==chr_list[n], c("position.start", "position.end"), drop=FALSE]
   message(paste(sep=" : ", chr_list[n], nrow(chr_loh_x)))
   chr_tumor_loh_y = c()
   if( nrow(chr_loh_x) !=0 ){
      for( i in 1:nrow(chr_loh_x) ){
         index = which( eLOH$position >= chr_loh_x$position.start[i] & eLOH$position <= chr_loh_x$position.end[i]
            & eLOH$chr == chr_list[n] )
         chr_loh_y = eLOH$tumor.baf[index] / eLOH$tumor.coverage[index]
            chr_loh_y_up = chr_loh_y[chr_loh_y>=0.5]
            chr_loh_y_down = chr_loh_y[chr_loh_y<0.5]

         chr_tumor_loh_y = rbind(chr_tumor_loh_y, c(mean(chr_loh_y_up), mean(chr_loh_y_down)))
      }
   }
  
   chr_loh_x = chr_loh_x + chr_length_cum[n]
   tumor_loh_x_y_chr = cbind(chr_loh_x, chr_tumor_loh_y)

}
# write.table(tumor_loh_x_y, "/WORK/GR/zhanglanying/work/WGS.breastCaner.3/LOH/test.txt", sep="\t", eol="\n", quote=FALSE, col.names=FALSE, row.names=FALSE)

opar = par(mar=c(5,4,4,2)+0.1)
png( filename=png_file , width = 2300, height = 400 )
par(mar=c(5,6,2.5,2))
### SNP_BAF
# baf_col = rainbow(7,0.5)
# plot( baf_x , tumor_baf_y , xlim = c(0,chr_length_cum[24]) , ylim =c(0,1.05),
  #  xlab = '' , ylab = 'BAF' , xaxt = 'n'  , xaxs = 'i' , yaxs = 'i' ,
  #  col = rainbow(7,0.8)[6] , pch = 16 , cex = 0.4 , bty='n', cex.lab = 2.5 , cex.axis = 2)

source("/WORK/GR/zhanglanying/program/scripts/getColors.R")
# cols = getColors(23,100)
cols = rep("gray40", 23)
names(cols) = chr_list

index = which(baf_x[,2]=="chr1")
plot( baf_x[index,1] , tumor_baf_y[index,1] , xlim = c(0,chr_length_cum[24]) , ylim =c(0,1.05),
   xlab = '' , ylab = 'BAF' , xaxt = 'n'  , xaxs = 'i' , yaxs = 'i' ,
   col = cols[1] , pch = 16 , cex = 0.4 , bty='n', cex.lab = 2.5 , cex.axis = 2)
for(chr in chr_list[2:23]){
   par(new=TRUE)
   index = which(baf_x[,2]==chr)
   plot( baf_x[index,1] , tumor_baf_y[index,1] , xlim = c(0,chr_length_cum[24]) , ylim =c(0,1.05),
      xlab = '' , ylab = '' , xaxt = 'n'  , yaxt='n', xaxs = 'i' , yaxs = 'i' ,
      col = cols[chr] , pch = 16 , cex = 0.4 , bty='n', cex.lab = 2.5 , cex.axis = 2)
}


### plot dotted line
abline( h = seq(0.1,1,by=0.2) , lty = 5 ,col = 'gray') 
### plot loh line
tumor_loh_num = nrow( tumor_loh_x_y )
message(paste(sep=": ","LOH",tumor_loh_num))
message(paste(sep=": ", "LOH_file", nrow(loh)))
for(i in 1:tumor_loh_num){
   # lines( loh_x[i,] , rep(tumor_loh_baf_y[i],2) , col = "forestgreen" , lwd=3 )
   lines( tumor_loh_x_y[i,c(1,2)] , rep(tumor_loh_x_y[i,3],2) , col = "darkorange" , lwd=3 )
   lines( tumor_loh_x_y[i,c(1,2)] , rep(tumor_loh_x_y[i,4],2) , col = "darkorange" , lwd=3 )
   # message(tumor_loh_baf_y[i,1])
}

### plot chr_line and add chr_main
abline( v = chr_length_cum[2:23] , col = 'gray44')
  
temp_chr =seq(1,22,by=1)
temp_chr = c(temp_chr, "X")
chr_main_x = c()
for(i in 1:23){
   temp = chr_length_cum[i] + chr_length[i]/2
   chr_main_x = c(chr_main_x, temp)
}
text( x = chr_main_x , y = rep(1.04,23) , labels = temp_chr , pos = 1 , cex=2)
dev.off()
par(opar)
