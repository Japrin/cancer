args = commandArgs(TRUE)
cnv_file = args[1]
gff_file = args[2]

cnv = read.table(cnv_file ,sep="\t", header=TRUE, as.is=TRUE)
index = which(cnv$copy.number==1 | cnv$copy.number==3)
cnv$copy.number[which(cnv$copy.number==1)] = rep("deletion", length(which(cnv$copy.number==1)))
cnv$copy.number[which(cnv$copy.number==3)] = rep("duplication", length(which(cnv$copy.number==3)))
out_9 = paste(sep="", 
   "copy.number=", cnv[index,8] , ";",
   "targeted.base=", cnv[index,5], ";",
   "coverage=", cnv[index,4], ";",
   "num.mark=", cnv[index,6], ";",
   "lower.cutoff=", cnv[index,9], ";", 
   "upper.cutoff=", cnv[index,10], ";",
   "log.ratio=", cnv[index,11], ";",
   "spec=", cnv[index,13], ";",
   "sens=", cnv[index,14])
out = cbind(cnv[index,1], rep("ExomeCNV",length(index)), cnv$copy.number[index], cnv[index,2], cnv[index,3],
   rep(".",length(index)), rep(".",length(index)), rep(".",length(index)) , out_9)
write.table(out, gff_file , quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE)
