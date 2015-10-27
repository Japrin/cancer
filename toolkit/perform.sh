#!/bin/bash

if [ $# -lt 1 ];then
	echo "usage perform <infile>"
	exit
fi

ifile=$1

#aln     aln     primary aln_NHDE0075_Run1.L8    15:02:58        54178   11911   8777
printf "ID1\tID2\tStage\tjob\twall_time\twall_time(second)\tCPU_time(second)\tMax Mem(MB)\n"
perl -e '@L=<>;for(my $i=0;$i<@L;$i++) { if($L[$i]=~/^Successful jobs/){$beg=$i}elsif($L[$i]=~/^Final summary/){$end=$i}} for(my $i=$beg+1;$i<$end;$i++){print $L[$i]} ' $ifile \
	| perl -ane 's/^somatic_//; if(/comparison/){next} /(.+?)_(.+?) \((.+?)\): (.+?) \((\d+) sec, (\d+)\/\d+ MB\)/; $ID2=(split /_/,$2)[0]; $ID1=(split /[AB]/,$ID2)[0];  printf "$ID1\t$ID2\t$1\t$2\t$4\t$5\t$6\n";' \
	| sort -k 3,3 \
	| perl -F"\t" -ane 'chomp @F;($h,$m,$s)=split /:/,$F[4];$wt=$h*60*60+$m*60+$s;splice @F,5,0,$wt;printf "%s\n",join("\t",@F);'
