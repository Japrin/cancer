### plot all_length for each chromosome
### ylimit
### Seg_num_markers
### CNV_highlight
#.libPaths(new="/WORK/GR/zhanglanying/01software/Rlib")
library(getopt)

spec = matrix( c( 
    "bin_file"    , "b" , 1 , "character" , "bin_log_ratio_file"         ,
    "seg_file"    , "s" , 1 , "character" , "segment_log_ratio_file"     , 
    "png_file"    , "p" , 1 , "character" , "png_file"                   ,

    "num_markers" , "n" , 2 , "integer"   , "CNV num_markers(default 8)" ,

    "y_bottom"    , "d" , 2 , "double"    , "y-axis bottom value"        ,
    "y_top"       , "u" , 2 , "double"    , "y-axis top value"           ,
  
    "cnv_file"    , "c" , 2, "character"  , "highlight CNV file"         ,

    "help"        , "h" , 0 , "logical"   , "help"
   ), ncol=5 , byrow=TRUE ) 

opt = getopt(spec , opt = commandArgs(TRUE))

if( !is.null(opt$help) ) {
   cat(getopt(spec , usage=TRUE))
   q(status=1)
}

if( is.null(opt$bin_file) || is.null(opt$seg_file) || is.null(opt$png_file) ){
   message("Error: 'bin_file' , 'seg_file' or 'png_file' is null")
   q()
}

if( is.null(opt$num_markers) ){ opt$num_markers = 8 }
if( is.null(opt$y_bottom)    ){ opt$y_bottom = -2   }
if( is.null(opt$y_top)       ){ opt$y_top = 2 }
if( opt$y_bottom >= opt$y_top ){
   message("Error: y_bottom >= y_top ")
   q()
}

# Bin_called_file_dir = "/PROJ/GR/HUMAN/zhujinde/cellline_zhujinde/all/cnv/NHDE0129/NHDE0129.somaticSNA.varscan.copynumber.called"
# Seg_called_file_dir = "/PROJ/GR/HUMAN/zhujinde/cellline_zhujinde/all/cnv/NHDE0129/NHDE0129.somaticSNA.varscan.copynumber.called.segment"
# png_file_name="/WORK/GR/zhanglanying/CNV/plot_CNV_22/NHDE0129.new.test.png"

### read Bin_called_data
Bin_called_data =  read.table( opt$bin_file , sep = "\t" , header = TRUE )
# colnames(Bin_called_data) = c( 'chrom' , 'chr_start' , 'chr_stop' , 'num_positions' , 
#   'normal_depth' , 'tumor_depth' , 'adjusted_log_ratio' , 'gc_content' , 'region_call' , 'raw_ratio' )

### read Seg_called_data
Seg_called_data = read.table( opt$seg_file , sep = "\t" , header = TRUE , as.is=TRUE )

### read cnv_file
if( !is.null(opt$cnv_file) ){
   cnv_called_data = read.table( opt$cnv_file , sep="\t" , header=TRUE , as.is=TRUE )
}

### chr base_information
chr_list = c(1:22,'X')
chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 
133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560 )
names(chr_length)=chr_list

x_limit = 0
Bin_logR_x_y = c()

Seg_logR_x = c()
Seg_logR_y = c()
Seg_num_markers = c()

cnv_logR_x = c()
cnv_logR_y = c()

chr_line_x = c()
chr_main_x = c()

for(chr in 1:length(chr_list)){
   chr_Bin_called_data = Bin_called_data[Bin_called_data$chrom==chr_list[chr],,drop=FALSE]
   chr_Seg_called_data = Seg_called_data[Seg_called_data$chrom==chr_list[chr],,drop=FALSE]   
   if( !is.null(opt$cnv_file) ){
      chr_cnv_called_data = cnv_called_data[cnv_called_data$chrom==chr_list[chr],,drop=FALSE]
   }
   
   if( nrow(chr_Bin_called_data) != 0 ){
      Bin_logR_x = chr_Bin_called_data[,'chr_start'] + x_limit
      Bin_logR_x_y = rbind( Bin_logR_x_y , cbind( Bin_logR_x , chr_Bin_called_data[,'adjusted_log_ratio'] ) )  
   }

   if( nrow(chr_Seg_called_data) != 0 ){
      Seg_logR_x = rbind( Seg_logR_x , (chr_Seg_called_data[,c('loc.start','loc.end'),drop=FALSE] + x_limit) )
      Seg_logR_y = rbind( Seg_logR_y , matrix(rep(chr_Seg_called_data[,'seg.mean'],2) , ncol=2) )
      Seg_num_markers = c( Seg_num_markers , chr_Seg_called_data[,'num.mark'])
   }

   if( !is.null(opt$cnv_file) ){
      if( nrow(chr_cnv_called_data) !=0 ){
         cnv_logR_x = rbind( cnv_logR_x , (chr_cnv_called_data[,c('chr_start','chr_stop'),drop=FALSE] + x_limit) )
         cnv_logR_y = rbind( cnv_logR_y , matrix(rep(chr_cnv_called_data[,'seg_mean'],2) , ncol=2) )
      }
   }
   
   chr_main_x = c( chr_main_x , mean(c(0,chr_length[chr]))+x_limit )
   x_limit = x_limit + chr_length[chr]
   chr_line_x = c( chr_line_x , x_limit )
}

### y_limit
if( is.null(opt$y_bottom) & is.null(opt$y_top) ){
   if( max(Seg_logR_y[,2])>=2 || min(Seg_logR_y[,2])<=(-2) ){
      opt$y_bottom = min( c(Bin_logR_x_y[,2] , Seg_logR_y[,2]) )
      opt$y_top = max( c(Bin_logR_x_y[,2] , Seg_logR_y[,2]) )
   }
}

opar = par(mar=c(5,4,4,2)+0.1)
png( filename = opt$png_file , width = 2300, height = 400 )
par(mar=c(1,6,1,1))
### plot Bin_logR point
plot( Bin_logR_x_y[,1] , Bin_logR_x_y[,2] , xlim = c(0,x_limit) , ylim = c(opt$y_bottom , opt$y_top*1.05) , 
   xlab = '' , ylab = 'Log Ratio' , cex.lab = 2.5 , cex.axis = 2 ,xaxt = 'n' , xaxs = 'i' , yaxs = 'i' ,
   col = 'gray90' , pch = 16 , cex = 0.4 , bty='n')
### plot dotted line
if( opt$y_bottom >= (-0.5) ){
   temp_1 = 0
}else{
   temp_1 = seq( 0 , opt$y_bottom + 0.5 , -0.5 )
}
if( opt$y_top <= 0.5 ){
   temp_2 = 0
}else{
   temp_2 = seq( 0 , opt$y_top - 0.5 , 0.5 )
}

dotted_line_y = c( temp_1 , temp_2 )
abline( h = dotted_line_y , lty = 5 ,col = 'gray')

### plot Seg_logR line
Seg_num = nrow( Seg_logR_x )

if( is.null(opt$cnv_file) ){
   for(i in 1:Seg_num){
      if( Seg_logR_y[i,1] >= 0.3 && Seg_num_markers[i] >= opt$num_markers ){
         Seg_line_col = 'red'
      }else{
         if( Seg_logR_y[i,1] <= -0.3 && Seg_num_markers[i] >= opt$num_markers ){
            Seg_line_col = 'forestgreen'
         }else{
            Seg_line_col = 'gold'
         }
      }
      lines( Seg_logR_x[i,] , Seg_logR_y[i,] , col = Seg_line_col , lwd=3 )
   }
}else{
   for(i in 1:Seg_num){
      lines( Seg_logR_x[i,] , Seg_logR_y[i,] , col = "gold" , lwd=2 )
   }

   cnv_num = nrow( cnv_logR_x )
   for(i in 1:cnv_num){
      if( cnv_logR_y[i,1] >= 0.3 ){
         cnv_line_col = 'red'
      }else{
         if( cnv_logR_y[i,1] <= -0.3){
            cnv_line_col = "forestgreen"
         }else{
            cnv_line_col = "gold"
         }
      }
      lines( cnv_logR_x[i,] , cnv_logR_y[i,] , col = cnv_line_col , lwd=3 )
   }
}
### plot chr_line and add chr_main
chr_num = length( chr_list )
chr_line_x=chr_line_x[-chr_num]
abline( v = chr_line_x , col = 'gray44')  
chr_main = chr_list 
text( x = chr_main_x , y = rep(opt$y_top*1.04,chr_num) , labels = chr_main , pos = 1  , cex=2 )
dev.off()
par(opar)

message(paste("Bin_log_ratio_file:   ",opt$bin_file,sep=""))
message(paste("Seg_log_ratio_file:   ",opt$seg_file,sep=""))
if( !is.null(opt$cnv_file) ){
   message(paste("CNV_log_ratio_file:   ",opt$cnv_file,sep=""))
}
message(paste("y_axis_range:   ",'[',opt$y_bottom,' , ',opt$y_top,']',sep=""))
if( is.null(opt$cnv_file) ){
   message(paste("Seg_num_markers:   ",opt$num_markers,sep=""))
}
message(paste("png_file:   ",opt$png_file,sep=""))
# message(cnv_num)
message("")

