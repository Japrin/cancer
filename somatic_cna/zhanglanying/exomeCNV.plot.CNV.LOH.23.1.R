### For ExomeCNV result
### plot all_length for each chromosome

args = commandArgs(TRUE)
if( length(args)==5 || length(args)==7 ){
   exon_lrr_txt = args[1]
   segment_lrr_txt = args[2]

   eLOH_txt = args[3]
   loh_txt = args[4]

   png_file = args[5]

   if(length(args)==7){
      if(as.numeric(args[6]) >=as.numeric(args[7])){
         message("Usage: Rscript plot_cnv_23.3.R <.exon.lrr.txt> <.segment.lrr.txt> <.eLOH.txt> <.loh.txt> <png_file> [y_floor_limit(default -2) y_ceiling_limit(default 2)]")
         stop()
      }
   }
}else{
   message("Usage: Rscript plot_cnv_23.3.R <.exon.lrr.txt> <.segment.lrr.txt> <.eLOH.txt> <.loh.txt> <png_file> [y_floor_limit(default -2) y_ceiling_limit(default 2)]")
   stop()
}

# exon_lrr_txt = "/PROJ/GR/HUMAN/intra/yongyong_cnv_20130927/analysis_0928/out_0930/YY15/somatic/YY15-1/TR_cnv/YY15.YY15-1.cnv.exon.lrr.txt"
# segment_lrr_txt = "/PROJ/GR/HUMAN/intra/yongyong_cnv_20130927/analysis_0928/out_0930/YY15/somatic/YY15-1/TR_cnv/YY15.YY15-1.cnv.segment.lrr.txt"
### chr base_information
chr_list = paste('chr',1:22,sep="")
chr_list = c(chr_list , 'chrX')
chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516,
   133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560 )
names(chr_length)=chr_list

if(file.info(segment_lrr_txt)$size !=0 ){
   ### read Exon_called_data
   Exon_called_data =  read.table( exon_lrr_txt , sep = "\t" , header = FALSE )
   colnames(Exon_called_data) = c( 'chrom' , 'chr_start' , 'chr_stop' , 'log_ratio' )

   ### read Seg_called_data
   Seg_called_data = read.table( segment_lrr_txt , sep = "\t" , header = FALSE )
   colnames(Seg_called_data) = c( 'chrom' , 'loc.start' , 'loc.end' , 'log_ratio' )

   x_limit = 0
   Exon_logR_x_y = c()
   Seg_logR_x = c()
   Seg_logR_y = c()
   chr_line_x = c()
   chr_main_x = c()

   for(chr in 1:length(chr_list)){
      chr_Exon_called_data = Exon_called_data[Exon_called_data$chrom==chr_list[chr],,drop=FALSE]
      chr_Seg_called_data = Seg_called_data[Seg_called_data$chrom==chr_list[chr],,drop=FALSE]

      if( nrow(chr_Exon_called_data) != 0 ){     
         Exon_logR_x = chr_Exon_called_data[,'chr_start'] + x_limit
         Exon_logR_x_y = rbind( Exon_logR_x_y , cbind( Exon_logR_x , chr_Exon_called_data[,'log_ratio'] , rep(chr_list[chr], length(Exon_logR_x)) ) )
      }

      if( nrow(chr_Seg_called_data) !=0 ){
         Seg_logR_x = rbind( Seg_logR_x , (chr_Seg_called_data[,c('loc.start','loc.end'),drop=FALSE] + x_limit) ) 
         Seg_logR_y = rbind( Seg_logR_y , matrix(rep(chr_Seg_called_data[,'log_ratio'],2) , ncol=2) ) 
      } 

      chr_main_x = c( chr_main_x , mean(c(0,chr_length[chr]))+x_limit )
      x_limit = x_limit + chr_length[chr]
      chr_line_x = c( chr_line_x , x_limit )
   }

}


if(file.info(loh_txt)$size !=0 ){
   #########  BAF  #########
   #########  BAF  #########
   eLOH =  read.table(eLOH_txt, sep="\t", header=TRUE, as.is=TRUE)
   chrY_index = which(eLOH$chr=="chrY")
   if(length(chrY_index)>0){
      eLOH = eLOH[-chrY_index,]
   }

   loh = read.table(loh_txt, sep="\t", header=TRUE, as.is=TRUE)
   loh = loh[(loh$LOH=="TRUE" & loh$num.mark>=4),]

   normal_baf_y = eLOH$normal.baf / eLOH$normal.coverage
   # tumor_baf_y = eLOH$tumor.baf / eLOH$tumor.coverage

   # normal_loh_baf_y = loh$normal.baf / loh$normal.coverage
   # tumor_loh_baf_y = loh$tumor.baf / loh$tumor.coverage 

   x_limit = 0
   baf_x = c()
   tumor_baf_y = c()

   tumor_loh_x = c()
   tumor_loh_y = c()

   chr_main_x = c()
   chr_line_x = c()

   sum = 0
   ### SNP_BAF

   for(chr in 1:length(chr_list)){
      chr_baf_x = eLOH[eLOH$chr==chr_list[chr],"position"]
      chr_baf_y = eLOH[eLOH$chr==chr_list[chr],"tumor.baf"] / eLOH[eLOH$chr==chr_list[chr],"tumor.coverage"]
      chr_loh_x = loh[loh$chr==chr_list[chr], c("position.start", "position.end"),drop=FALSE]
   
      if( length(chr_baf_x) != 0 ){
         baf_x = rbind( baf_x, cbind(chr_baf_x + x_limit, rep(chr_list[chr], length(chr_baf_x)) ) )
         tumor_baf_y = rbind(tumor_baf_y, cbind(chr_baf_y, rep(chr_list[chr], length(chr_baf_y)) ))
      }

      if( nrow(chr_loh_x) !=0 ){ 
         for( i in 1:nrow(chr_loh_x) ){
            index = which( eLOH$position >= chr_loh_x$position.start[i] & eLOH$position <= chr_loh_x$position.end[i]
               & eLOH$chr == chr_list[chr] )
            sum = length(index) + sum
            chr_loh_y = eLOH$tumor.baf[index] / eLOH$tumor.coverage[index]
         
               chr_loh_y_up = chr_loh_y[chr_loh_y>=0.5]
               chr_loh_y_down = chr_loh_y[chr_loh_y<=0.5]
        
            tumor_loh_y = rbind(tumor_loh_y, c(mean(chr_loh_y_up), mean(chr_loh_y_down)))
         }
         tumor_loh_x = rbind( tumor_loh_x , chr_loh_x + x_limit)
      }
      chr_main_x = c( chr_main_x , mean(c(0,chr_length[chr]))+x_limit )
      x_limit = x_limit + chr_length[chr]
      chr_line_x = c( chr_line_x , x_limit )
   }
   
   message(sum)
   message(sum(loh$num.mark))

}


if(file.info(segment_lrr_txt)$size == 0 & file.info(loh_txt)$size ==0){
   png( filename=png_file , width = 2300, height = 800 )
   dev.off()
}

if(file.info(segment_lrr_txt)$size != 0 | file.info(loh_txt)$size !=0){

   opar = par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
   if(file.info(segment_lrr_txt)$size != 0 & file.info(loh_txt)$size !=0){
      png( filename=png_file , width = 2300, height = 800 )
      par(mfrow=c(2,1))
   }else{
      png( filename=png_file , width = 2300, height = 400 )
   }
   #########  CNV  #########
   #########  CNV  #########
   if(file.info(segment_lrr_txt)$size != 0){
      # cols = rep("gray90", 23)
      cols = rep("black", 23)
      if(file.info(loh_txt)$size !=0){
         par(mar=c(0,6,4,2))
      }else{
         par(mar=c(2,6,2,2))
      }

      # y_limit = range( c(Exon_logR_x_y[,2],Seg_logR_y[,2]) ) 
      if(length(args)==7){
         y_limit = as.numeric(c(args[6],args[7]))
      }else{
         if( max(Seg_logR_y[,2])>=2 | min(Seg_logR_y[,2])<=(-2) ){
            y_limit = range( as.numeric(Exon_logR_x_y[,2]) )
         }else{
            y_limit = c(-2,2)
         }
      }
      # col = cols[1]
      plot( Exon_logR_x_y[,1] , Exon_logR_x_y[,2] , xlim = c(0,x_limit) , ylim = c(y_limit[1],y_limit[2]*1.05), 
         xlab = '' , ylab = 'Log Ratio' , xaxt = 'n' , xaxs = 'i' , yaxs = 'i' ,
         col = cols[1], pch = 16 , cex = 0.4 , bty='n', cex.lab = 2.5 , cex.axis = 2)
      
      ### plot dotted line
      if( y_limit[1] >= (-0.5) ){
         temp_1 = 0
      }else{
         temp_1 = seq( 0 , y_limit[1]+0.5 , -0.5 )
      }

      if( y_limit[2] <= 0.5 ){
         temp_2 = 0
      }else{
         temp_2 = seq( 0 , y_limit[2]-0.5 , 0.5 )
      }

      dotted_line_y = c( temp_1 , temp_2 )
      abline( h = dotted_line_y , lty = 5 ,col = 'gray')

      ### plot Seg_logR line
      Seg_num = nrow( Seg_logR_x )
      for(i in 1:Seg_num){
         if( Seg_logR_y[i,1] >= 0.3 ){
            Seg_line_col = 'red'
            # Seg_line_col = "darkorange"
            # Seg_line_col = "red"
         }else{
            if( Seg_logR_y[i,1] <= -0.3 ){
               # Seg_line_col = 'forestgreen'
               # Seg_line_col = "green"
               # Seg_line_col = "darkorange"
               Seg_line_col = "red"
            }else{
               # Seg_line_col = 'darkorange'
               Seg_line_col = "red"
            }
         }
         lines( Seg_logR_x[i,] , Seg_logR_y[i,] , col = Seg_line_col , lwd=3)
      }
      ### plot chr_line and add chr_main
      chr_num = length( chr_list )
      chr_line_x=chr_line_x[-chr_num]
      abline( v = chr_line_x , col = 'gray44')  
      temp_chr =seq(1,22,by=1)
      temp_chr = c(temp_chr, "X")
      text( x = chr_main_x , y = rep(y_limit[2]*1.04,23) , labels = temp_chr , pos = 1 , cex=2)
      # abline( h = 0 , lty = 5 )
      message("CNV finished")
   }

   if(file.info(loh_txt)$size != 0){
      ### plot LOH ###
      ### plot LOH ###
      if(file.info(segment_lrr_txt)$size !=0){
         par(mar=c(5,6,2.5,2))
      }else{
         par(mar=c(2,6,2,2))
      }
      #  cols = rep("gray80", 23)
      cols = rep("black", 23)

      if(file.info(segment_lrr_txt)$size == 0){
         hah = 1.1
      }else{
         hah = 1
      }
      # col = cols[1]
      plot( baf_x[,1] , tumor_baf_y[,1] , xlim = c(0,x_limit) , ylim =c(0,hah),
         xlab = '' , ylab = 'BAF' , xaxt = 'n'  , xaxs = 'i' , yaxs = 'i' ,
         col = cols[1] , pch = 16 , cex = 0.4 , bty='n', cex.lab = 2.5 , cex.axis = 2)

      ### plot dotted line
      abline( h = seq(0.1,1,by=0.2) , lty = 5 ,col = 'gray') 
      ### plot loh line
      if( length(tumor_loh_x)>0 ){
         tumor_loh_num = nrow( tumor_loh_x )
         for(i in 1:tumor_loh_num){
            # lines( loh_x[i,] , rep(tumor_loh_baf_y[i],2) , col = "forestgreen" , lwd=3 )
            # darkorange1
            # lines( tumor_loh_x[i,] , rep(tumor_loh_y[i,1],2) , col = "red" , lwd=3 )
           #  lines( tumor_loh_x[i,] , rep(tumor_loh_y[i,2],2) , col = "red" , lwd=3 )
            # message(tumor_loh_baf_y[i,1])
         }
      }
      ### plot chr_line and add chr_main
      abline( v = chr_line_x , col = 'gray44')
      if(file.info(segment_lrr_txt)$size == 0){
         temp_chr =seq(1,22,by=1)
         temp_chr = c(temp_chr, "X")
         text( x = chr_main_x , y = rep(hah,23) , labels = temp_chr , pos = 1 , cex=2)
      }
   }
   dev.off()
   par(opar)
}

if(file.info(segment_lrr_txt)$size == 0 & file.info(loh_txt)$size ==0){
   png( filename=png_file , width = 400, height = 400 )
   plot(1:10,1:10,type="n",xaxt='n',yaxt='n',xlab="",ylab="",bty='n')
   dev.off()
}
