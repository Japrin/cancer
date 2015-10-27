args = commandArgs(T)
cnv_cbs_txt = args[1]
out_file = args[2]

da = read.table(cnv_cbs_txt, sep="\t", header=T, as.is=T)

mark = c(2,8,20)
n = length(mark) + 1
out = c("del_number", "del_length_of_target", "del_length_of_all", 
        "dup_number", "dup_length_of_target", "dup_length_of_all")
for(i in 1:n){
   if(i==1){
      index_del = which(da$num.mark>0 & da$num.mark<=mark[i] & da$copy.number==1)
      index_dup = which(da$num.mark>0 & da$num.mark<=mark[i] & da$copy.number==3)
   }else{
      if(i==n){
         index_del = which(da$num.mark>mark[i-1] & da$copy.number==1)
         index_dup = which(da$num.mark>mark[i-1] & da$copy.number==3)    
      }else{
         index_del = which(da$num.mark>mark[i-1] & da$num.mark<=mark[i] & da$copy.number==1)
         index_dup = which(da$num.mark>mark[i-1] & da$num.mark<=mark[i] & da$copy.number==3)
      }  
   }

   num_del = length(index_del)
   num_dup = length(index_dup)
   
   if(length(num_del)>0){
      sum_all_del = sum( da$probe_end[index_del] - da$probe_start[index_del] + 1 )
      sum_target_del = sum( da$targeted.base[index_del] )
   }else{
      sum_all_del = 0
      sum_target_del = 0
   }

   if(length(num_dup)>0){
      sum_all_dup = sum( da$probe_end[index_dup] - da$probe_start[index_dup] + 1)
      sum_target_dup = sum( da$targeted.base[index_dup] )
   }else{
      sum_all_dup = 0
      sum_target_dup = 0
   }

   temp = c(num_del , sum_target_del, sum_all_del, 
            num_dup , sum_target_dup, sum_all_dup)

   out = cbind(out, temp)
}

label = paste(sep="", 'num.mark<=',mark)
temp = c(0, mark[1:(n-2)] )
label = paste(sep="", temp,'<', label)
label = c(label , paste(sep="", mark[n-1],"<num.mark"))
label = c("category", label)
colnames(out) = label

write.table(out , out_file , sep="\t", eol="\n" , quote=F, row.names=F, col.names=T)
