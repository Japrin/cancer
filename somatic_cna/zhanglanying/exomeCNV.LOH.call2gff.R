args = commandArgs(TRUE)
loh_file = args[1]
gff_file = args[2]

loh = read.table(loh_file ,sep="\t", header=TRUE, as.is=TRUE)
index = which(loh$LOH=="TRUE")
index_2 = which(loh$LOH=="TRUE" & loh$num.mark==1)
loh[index_2,2] = loh[index_2,2] - 1 

out_9 = paste(sep="", 
   "num.mark=", loh[index,4], ";",
   "normal.coverage=", loh[index,5], ";",
   "normal.baf=", loh[index,6], ";",
   "tumor.coverage=", loh[index,7], ";",
   "tumor.baf=", loh[index,8], ";",
   "p_value=", loh[index,10])
out = cbind(loh[index,1], rep("ExomeCNV",length(index)), rep("LOH", length(index)), loh[index,2], loh[index,3],
   rep(".",length(index)), rep(".",length(index)), rep(".",length(index)) , out_9)
write.table(out, gff_file , quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE)
