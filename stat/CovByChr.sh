#!/bin/bash

if [ $# -lt 2 ]
then
	echo "Usage: $0 <in> <out>"
	exit 1
fi

ifile=$1
ofile=$2

tDir=.
tfile=$tDir/tmp.123456.ooo

#sed '1,1d' $ifile | perl -F"\t" -ane '{chomp; chomp @F; @a=split /[:-]/,$F[0];if(@a<3){push @a,$a[1];} $F[0]=join("\t",@a);$len=$a[2]-$a[1]+1;print "$F[0]\t$len\t$F[3]\t$F[4]\t$F[8]\n"; }' > $tfile
cut -f 1,4,5,9 $ifile | sed -e '1,1d' -e 's/[:-]/\t/g' | awk -v OFS="\t" '{$3=$3"\t"($3-$2+1);print $0}' > $tfile

printf "chr\tlength\ttotal_coverage\tmean_coverage\tcovered_bases\tprop_covered_bases\n" > $ofile
for i in {1..22} X Y
do
	awk -v OFS="\t" ' $1=="'$i'"{chr=$1;l+=$4;if($5>0 && $6>0) { l2+=$5/$6;b+=$5;c+=int($4*$7/100+0.5)} else{l2+=$4} }END{print chr,l,b,b/l2,c,c/l}' $tfile 
done >> $ofile
rm $tfile
