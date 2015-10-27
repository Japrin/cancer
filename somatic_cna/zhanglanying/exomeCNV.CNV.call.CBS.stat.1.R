args = commandArgs(T)
cnv_cbs_txt = args[1]
out_file = args[2]

da = read.table(cnv_cbs_txt, sep="\t", header=T, as.is=T)

out = c("del_number", "del_length_of_target", "del_length_of_all", 
        "dup_number", "dup_length_of_target", "dup_length_of_all")

index_del = which(da$num.mark>=2 & da$copy.number=="deletion")
index_dup = which(da$num.mark>=2 & da$copy.number=="duplication")

   num_del = length(index_del)
   num_dup = length(index_dup)
   

# Start   End
   if(length(num_del)>0){
      sum_all_del = sum( da$End[index_del] - da$Start[index_del] + 1 )
      sum_target_del = sum( da$targeted.base[index_del] )
   }else{
      sum_all_del = 0
      sum_target_del = 0
   }

   if(length(num_dup)>0){
      sum_all_dup = sum( da$End[index_dup] - da$Start[index_dup] + 1)
      sum_target_dup = sum( da$targeted.base[index_dup] )
   }else{
      sum_all_dup = 0
      sum_target_dup = 0
   }

   temp = c(num_del , sum_target_del, sum_all_del, 
            num_dup , sum_target_dup, sum_all_dup)

   out = cbind(out, temp)

write.table(out , out_file , sep="\t", eol="\n" , quote=F, row.names=F, col.names=F)
