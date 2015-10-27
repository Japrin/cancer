### Rscript exomeCNV.LOH.R /PROJ/GR/HUMAN/sci/WES.shiguanai.12/02.analysis/EC952/EC952.M/var/EC952.EC952.M.sample_GATK.snp.call.reformated.vcf.gz /PROJ/GR/HUMAN/sci/WES.shiguanai.12/02.analysis/EC952/EC952.Dan/var/EC952.EC952.Dan.sample_GATK.snp.call.reformated.vcf.gz EC952.Dan /WORK/GR/zhanglanying/program/analysis-pipeline/12.ExomeCNV/test

library(ExomeCNV,lib.loc="/PROJ/CR/share/pipeline/cancerPipeline/Bin/ExomeCNV.v1/Rlib")

args = commandArgs(TRUE)
normal_vcf = args[1]
tumor_vcf = args[2]
tumor_sample_id = args[3]
output_dir = args[4]


a=0
if(a==1){
   out_merge_vcf = paste(sep="", output_dir,'/',tumor_sample_id,'.merge.vcf' )

   system( paste("/usr/java/latest/bin/java -jar /PROJ/GR/share/Software/medinfo/01bin/gatk/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar -T CombineVariants -R /PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/human_g1k_v37_decoy.fasta --variant ", normal_vcf, " --variant ", tumor_vcf, " -o ", out_merge_vcf, " -genotypeMergeOptions UNIQUIFY", sep="") )

   out_normal_baf = paste(sep='', output_dir, '/', tumor_sample_id, '.normal.baf')
   out_tumor_baf = paste(sep='',output_dir, '/', tumor_sample_id, '.tumor.baf')

   temp = read.table(out_merge_vcf, nrows=300, sep="*", comment.char="*",as.is=TRUE)
   index = grep("#",temp[,1])
   index = max(index)
   temp = temp[index,1]

   temp = unlist(strsplit(temp,split="\t"))
   tumor_index = grep("variant2",temp)
   message(tumor_index)
   if( tumor_index == 10){
      normal_col = 2
      tumor_col = 1
   }
   if( tumor_index == 11){
      tumor_col = 2
      normal_col = 1
   }
   message(paste(sep="", 'normal_col: ', normal_col))
   message(paste(sep="", 'tumor_col: ', tumor_col))
}

system( paste(sep="", "/PROJ/CR/zhanglanying/bin/ExomeCNV/ExomeCNV.human/exomeCNV.LOH.BAF.pl ", normal_vcf, " ", tumor_vcf, " ", tumor_sample_id, " ", output_dir) )
out_normal_baf = paste(sep='', output_dir, '/', tumor_sample_id, '.normal.baf')
out_tumor_baf = paste(sep='',output_dir, '/', tumor_sample_id, '.tumor.baf')
message(paste(sep="","normal_BAF: ",out_normal_baf))
message(paste(sep="","tumor_BAF: ",out_tumor_baf))

### LOH_analysis
normal_baf = out_normal_baf
tumor_baf = out_tumor_baf
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
the.strata=make.loh.strata(normal,tumor)
the.strata.ls=list(the.strata)

eLOH = LOH.analyze(normal, tumor, strata=the.strata.ls, alpha=0.01,  method="deviation.half.norm")

set.seed(50)
the.loh = multi.LOH.analyze(normal, tumor, all.loh.ls=list(eLOH), test.alpha=0.001, method="deviation.half.norm", sdundo=c(1,2), alpha=c(0.05,0.01))

loh_png = paste(sep="",loh_output_prefix, '.LOH.png')
png(loh_png, res=70, width=2000, height=1200, pointsize=16)
do.plot.loh(the.loh, normal, tumor, "two.sample.fisher", plot.style="dev", num.mark.min=3)
dev.off()
message(paste(sep="",'LOH_png: ',loh_png))

loh_out = paste(sep="", loh_output_prefix, '.LOH.txt')
write.table(the.loh, loh_out, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
message(paste(sep="",'LOH_txt: ', loh_out))

eLOH_out = paste(sep="",output_dir, '/', tumor_sample_id, '.eLOH.txt')
write.table(eLOH, eLOH_out, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Admixture Rate
index = which(the.loh$num.mark>=3 & the.loh$LOH=="TRUE")
seg_c = 1 - 2*mean( abs(the.loh$tumor.baf[index]/the.loh$tumor.coverage[index] - 0.5) )

message(sum(the.loh$num.mark[index]))
bin_index = c()
for(i in index){
   temp = which( eLOH$position >= the.loh$position.start[i] & eLOH$position <= the.loh$position.end[i] & eLOH$chr == the.loh$chr[i] )
   bin_index = c(bin_index, temp)
}
message(length(bin_index))
seg_bin_c = 1 - 2*mean( abs(eLOH$tumor.baf[bin_index]/eLOH$tumor.coverage[bin_index] -0.5) )

index = which(eLOH$LOH=="TRUE")
bin_c_1 = 1 - 2*mean( abs(eLOH$tumor.baf[index]/eLOH$tumor.coverage[index] -0.5) )

eLOH = LOH.analyze(normal, tumor, strata=the.strata.ls, alpha=0.001, method="deviation.half.norm")
index = which(eLOH$LOH=="TRUE")
bin_c_2 = 1 - 2*mean( abs(eLOH$tumor.baf[index]/eLOH$tumor.coverage[index] -0.5) )

eLOH = LOH.analyze(normal, tumor, strata=the.strata.ls, alpha=0.0001, method="deviation.half.norm")
index = which(eLOH$LOH=="TRUE")
bin_c_3 = 1 - 2*mean( abs(eLOH$tumor.baf[index]/eLOH$tumor.coverage[index] -0.5) )

c_out = c("LOH$num.mark>=3(p<0.001)", seg_c , "LOH$num.mark>=3(p<0.001) & bin(p<0.01)", seg_bin_c,  "eLOH(p<0.01)", bin_c_1, "eLOH(p<0.001)", bin_c_2, "eLOH(p<0.0001)", bin_c_3)

c_out = as.matrix(c_out)
c_file = paste(sep="", output_dir, '/', tumor_sample_id, ".admixture_rate")

write.table(c_out, c_file, eol="\n", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
message(paste("Admixture Rate: ",c_file,sep=""))
