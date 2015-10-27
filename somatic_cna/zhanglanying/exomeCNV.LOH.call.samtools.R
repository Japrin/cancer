### Rscript exomeCNV.LOH.call.samtools.R /PROJ/GR/HUMAN/sci/WGS.breastCaner.3/zhengliangtao/YMY/snp/2711B/2711B.all.snp.vcf.gz /PROJ/GR/HUMAN/sci/WGS.breastCaner.3/zhengliangtao/YMY/snp/109PT-2/109PT-2.all.snp.vcf.gz 109PT /WORK/GR/zhanglanying/program/analysis-pipeline/12.ExomeCNV/test

### chr1    10109   .       A       T       40      .       Func=intergenic;Gene=NONE(dist=NONE),DDX11L1(dist=1765);SegDup=0.99;DP=117;VDB=3.924378e-03;RPB=1.590 163e-01;AF1=0.5;AC1=1;DP4=28,35,19,21;MQ=21;FQ=43;PV4=0.84,0.017,0.0031,1       GT:PL:DP:SP:GQ  0/1:70,0,207:103:1:73

.libPaths(new=c("/WORK/GR/zhanglanying/01software/Rlib"))
suppressMessages(library(doMC))
library(ExomeCNV,lib.loc="/WORK/GR/zhanglanying/01software/Rlib/test")

args = commandArgs(TRUE)
normal_vcf = args[1]
tumor_vcf = args[2]
tumor_sample_id = args[3]
output_dir = args[4]


# filter = args[5]
# message(filter)
# if(filter=="Y"){
#   system( paste(sep="", "/WORK/GR/zhanglanying/program/analysis-pipeline/ExomeCNV/exomeCNV.LOH.BAF.samtools.filter.pl ", normal_vcf, " ", tumor_vcf, " ", tumor_sample_id, " ", output_dir) )
# }
# if(filter=="N"){
#   system( paste(sep="", "/WORK/GR/zhanglanying/program/analysis-pipeline/ExomeCNV/exomeCNV.LOH.BAF.samtools.pl ", normal_vcf, " ", tumor_vcf, " ", tumor_sample_id, " ", output_dir) )
# }

system( paste(sep="", "/WORK/GR/zhanglanying/program/analysis-pipeline/ExomeCNV/exomeCNV.LOH.BAF.samtools.pl ", normal_vcf, " ", tumor_vcf, " ", tumor_sample_id, " ", output_dir) )

### LOH_analysis
normal_baf = paste(sep="",output_dir,'/',tumor_sample_id,".normal.baf")
tumor_baf = paste(sep="", output_dir, "/", tumor_sample_id, ".tumor.baf")
loh_output_prefix = paste(sep="", output_dir, '/' , tumor_sample_id)

normal = read.delim(normal_baf,header=TRUE)
tumor = read.delim(tumor_baf,header=TRUE)
del_index = c( which(is.na(normal$baf)=="TRUE") , which(is.na(tumor$baf)=="TRUE") , which(normal$coverage==0), which(tumor$coverage==0))
del_index = unique(del_index)
if(length(del_index)>1){
   normal = normal[-del_index,]
   tumor = tumor[-del_index,]
}
normal$chr = paste("chr",normal$chr,sep="")
tumor$chr = paste("chr",tumor$chr,sep="")

# eLOH = LOH.analyze(normal, tumor, alpha=0.05, method="two.sample.fisher")
# the.loh = multi.LOH.analyze(normal, tumor, all.loh.ls=list(eLOH), test.alpha=0.001, method="deviation.half.norm", sdundo=c(1,2), alpha=c(0.05,0.01))

###    normal   tumor
### 1. no_LOH   no_LOH
### 2. no_LOH      LOH
### 3.    LOH   no_LOH
### 4.    LOH      LOH
a=0
if(a==1){
  ### filter position which is 3 and 4.
  registerDoMC()
  chr.list = paste(sep="","chr",1:22)
  chr.list = c(chr.list,"chrX","chrY")
  temp_eLOH<- foreach(i=1:length(chr.list) , .combine=rbind , .inorder=TRUE) %dopar% {
     idx = which(normal$chr==chr.list[i])
     message(chr.list[i])
     # eloh = LOH.analyze(normal=normal[idx,], tumor=tumor[idx,], alpha=0.05, method="two.sample.fisher")
     eloh = LOH.analyze(normal=normal[idx,], tumor=tumor[idx,], alpha=0.01, method="only.normal")
  }
  message(paste(sep="", "total SNP count: ", nrow(temp_eLOH)))
  normal_LOH_index = which(temp_eLOH$LOH=="TRUE")
  temp_eLOH = temp_eLOH[-normal_LOH_index,,drop=FALSE]
  message(paste(sep="", "deletion position which is LOH in noraml sample: ", length(normal_LOH_index)))
  message(paste(sep="", "the count of SNP used in LOH analysis: ", nrow(temp_eLOH)))

  normal = temp_eLOH[,c(1,2,5,4)]
  colnames(normal) = c("chr", "position", "coverage", "baf")
  tumor = temp_eLOH[,c(1,2,7,6)]
  colnames(tumor) = c("chr", "position", "coverage", "baf")
}

### eLOH
### deviation.t
### deviation.half.norm
chr.list = paste(sep="","chr",1:22)
chr.list = c(chr.list,"chrX","chrY")
registerDoMC()
eLOH<- foreach(i=1:length(chr.list) , .combine=rbind , .inorder=TRUE) %dopar% {
   idx = which(normal$chr==chr.list[i])
   message(chr.list[i])

   the.strata=make.loh.strata(normal[idx,],tumor[idx,])
   the.strata.ls=list(the.strata)

   eloh = LOH.analyze(normal=normal[idx,], tumor=tumor[idx,], the.strata.ls, alpha=0.01, method="deviation.half.norm")
}

### deviation.half.norm
message("prepare CBS")
set.seed(50)
the.loh<- foreach(i=1:length(chr.list) , .combine=rbind , .inorder=TRUE) %dopar% {
   idx = which(normal$chr==chr.list[i])
   loh = multi.LOH.analyze(normal=normal[idx,], tumor=tumor[idx,], all.loh.ls=list(eLOH[which(eLOH$chr==chr.list[i]),,drop=FALSE]), 
      test.alpha=0.005, method="deviation.half.norm", sdundo=c(0,0), alpha=c(0.05,0.01) )
}

loh_png = paste(sep="",loh_output_prefix, '.LOH.png')
png(loh_png, res=70, width=2000, height=1200, pointsize=16)
do.plot.loh(the.loh, normal, tumor, "deviation.half.norm", plot.style="dev", num.mark.min=3)
dev.off()
message(paste(sep="",'LOH_png: ',loh_png))

loh_out = paste(sep="", loh_output_prefix, '.LOH.txt')
write.table(the.loh, loh_out, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
message(paste(sep="",'LOH_txt: ', loh_out))

eLOH_out = paste(sep="",output_dir, '/', tumor_sample_id, '.eLOH.txt')
write.table(eLOH, eLOH_out, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

message("prepare admixture rate")
# Admixture Rate
## 1/5
index = which(the.loh$num.mark>=3 & the.loh$LOH=="TRUE")
seg_c = 1 - 2*mean( abs(the.loh$tumor.baf[index]/the.loh$tumor.coverage[index] - 0.5) )
message("1/5")
message(sum(the.loh$num.mark[index]))

## 2/5
bin_index<- foreach(i=1:length(chr.list) , .combine=c , .inorder=TRUE) %dopar% {
   index = which(the.loh$chr==chr.list[i] & the.loh$num.mark>=3 & the.loh$LOH=="TRUE")
   message(paste(chr.list[i], ": ", length(index),sep=""))
   bin_index_chr = c()
   for(i in index){
      temp = which( eLOH$position >= the.loh$position.start[i] & eLOH$position <= the.loh$position.end[i] & eLOH$chr == the.loh$chr[i] )
      bin_index_chr = c(bin_index_chr, temp)
   }
   bin_index_chr = bin_index_chr
}
seg_bin_c = 1 - 2*mean( abs(eLOH$tumor.baf[bin_index]/eLOH$tumor.coverage[bin_index] -0.5) )
message("2/5")
message(length(bin_index))

## 3/5
index = which(eLOH$LOH=="TRUE")
bin_c_1 = 1 - 2*mean( abs(eLOH$tumor.baf[index]/eLOH$tumor.coverage[index] -0.5) )
message("3/5")

## 4/5
eLOH<- foreach(i=1:length(chr.list) , .combine=rbind , .inorder=TRUE) %dopar% {
   idx = which(normal$chr==chr.list[i])
   message(chr.list[i])
   the.strata=make.loh.strata(normal[idx,],tumor[idx,])
   the.strata.ls=list(the.strata)

   eloh = LOH.analyze(normal=normal[idx,], tumor=tumor[idx,], the.strata.ls, alpha=0.001, method="deviation.half.norm")
}
# eLOH = LOH.analyze(normal, tumor, alpha=0.001, method="two.sample.fisher")
index = which(eLOH$LOH=="TRUE")
bin_c_2 = 1 - 2*mean( abs(eLOH$tumor.baf[index]/eLOH$tumor.coverage[index] -0.5) )
message("4/5")

## 5/5
# eLOH = LOH.analyze(normal, tumor, alpha=0.0001, method="two.sample.fisher")
eLOH<- foreach(i=1:length(chr.list) , .combine=rbind , .inorder=TRUE) %dopar% {
   idx = which(normal$chr==chr.list[i])
   the.strata=make.loh.strata(normal[idx,],tumor[idx,])
   the.strata.ls=list(the.strata)
   eloh = LOH.analyze(normal=normal[idx,], tumor=tumor[idx,],the.strata.ls, alpha=0.0001, method="deviation.half.norm")
}
index = which(eLOH$LOH=="TRUE")
bin_c_3 = 1 - 2*mean( abs(eLOH$tumor.baf[index]/eLOH$tumor.coverage[index] -0.5) )
message("5/5")

### output_result
c_out = c("LOH$num.mark>=3(p<0.005)", seg_c , "LOH$num.mark>=3(p<0.005) & bin", seg_bin_c,  "eLOH(p<0.01)", bin_c_1, "eLOH(p<0.001)", bin_c_2, "eLOH(p<0.0001)", bin_c_3)

c_out = as.matrix(c_out)
c_file = paste(sep="", output_dir, '/', tumor_sample_id, ".admixture_rate")

write.table(c_out, c_file, eol="\n", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
message(paste("Admixture Rate: ",c_file,sep=""))
