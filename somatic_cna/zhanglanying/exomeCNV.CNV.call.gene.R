library(DNAcopy,lib.loc="/PUBLIC/software/public/R/v3.0.3/lib64/R/library")
library(ExomeCNV,lib.loc="/PROJ/CR/share/pipeline/cancerPipeline/Bin/ExomeCNV.v1/Rlib/")

# eCNV_file = paste(sep="",output_dir,'/',tumor_ID,'.eCNV.txt')
args = commandArgs(TRUE)
eCNV_file = args[1]
tumor_ID = args[2]
output_dir = args[3]

###########  GENE CNV
# eCNV = read.table("/WORK/GR/zhanglanying/program/analysis-pipeline/ExomeCNV/gene/EC949.Ca.exon.txt.gene", sep="\t", header=TRUE, as.is=TRUE)
system( paste(sep="", "perl /PROJ/CR/zhanglanying/bin/ExomeCNV/ExomeCNV.human/exomeCNV.CNV.eCNV2gene.pl ", eCNV_file) )
eCNV = read.table(paste(eCNV_file,'.gene',sep=""), sep="\t", header=TRUE, as.is=TRUE)

gene_region = read.table("/PROJ/CR/zhanglanying/bin/ExomeCNV/ExomeCNV.human/hg19_refGene.txt.gene", sep="\t", header=FALSE, as.is=TRUE, colClasses="character")

gene = unique(eCNV[,8])
if(length(which(is.na(gene))) > 0){
   gene = gene[-which(is.na(gene))]
}
message(paste(sep="","Gene: ", length(gene)))

out = c()
for(i in 1:length(gene)){
#   message(i)
   index = which(eCNV[,8]==gene[i])
   targeted_base = sum(eCNV$target.base[index])
   normal_average = sum(eCNV$normal.coverage[index]) / targeted_base
   tumor_average = sum(eCNV$tumor.coverage[index]) / targeted_base
   exon_num = length(index)  

   copy_0 = which(eCNV[index,"copy.number"]=="0")
   copy_1 = which(eCNV[index,"copy.number"]=="1")
   copy_2 = which(eCNV[index,"copy.number"]=="2")
   copy_3 = which(eCNV[index,"copy.number"]=="3")
   copy_NA = which(is.na(eCNV[index,"copy.number"]))
   copy_all = c( length(copy_0), length(copy_1), length(copy_2), length(copy_3), length(copy_NA)) 
   copy_index = c(0, 1, 2, 3, "NA")
   rank_copy = rank( copy_all )

   copy_max = which(rank_copy == max(rank_copy))

   if( (length(copy_1)>0 & length(copy_3)>0) || length(copy_max)>1 ){
      copy ="NA"
   }else{
      copy = copy_index[copy_max]
   }
   
   temp = c(as.matrix(gene_region[which(gene_region[,1]==gene[i]),])[1,], exon_num, targeted_base, normal_average, tumor_average, copy, copy_all)
   out = rbind(out, temp)
   eCNV = eCNV[-index,,drop=FALSE]
}

total.cov.normal = sum(as.numeric(out[,8])*as.numeric(out[,7]))
total.cov.tumor = sum(as.numeric(out[,9])*as.numeric(out[,7]))
logR = log2(as.numeric(out[,9])/as.numeric(out[,8])) + log2(total.cov.normal/total.cov.tumor)

normal.chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
logR = normalize.logR(logR, chr.median, out[,2] %in% normal.chrs)
out = cbind(out, logR)

colnames(out) = c("Gene", "Chr", "Transcription.start", "Transcription.end", "Exon.count", "Bin.count", 
   "Targeted.base", "Normal.average.coverage", "Tumor.average.coverage", "Copy.number",
   "Copy_0", "Copy_1", "Copy_2", "Copy_3", "Copy_NA","Log.ratio")
gene_CNV_file = paste(sep="" , output_dir , "/" , tumor_ID, ".CNV.gene.txt")
write.table(out, gene_CNV_file, eol="\n", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

index = which(out[,10]==1 | out[,10]==3)
out = cbind(out[index,2],out[index,3],out[index,4],out[index,16])
gene_CNV_file = paste(sep="" , output_dir , "/" , tumor_ID, ".CNV.gene.lrr.txt")
write.table(out, gene_CNV_file, eol="\n", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
