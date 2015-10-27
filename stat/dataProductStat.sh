#!/bin/bash

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <in> <out>"
	printf "\t<in> format: sample\tlane\t.stat file\n"
	exit 1
fi
ifile=$1
ofile=$2

while read s l f
do 
	#awk -F"\t" -f /PROJ/GR/share/medinfo.02pipeline/cancer/toolkit/transpose.awk $f \
	cat $f | transpose.sh \
		| awk -v OFS="\t" 'NR==1{print "Sample\tLane",$0}NR>1{print "'$s'""\t""'$l'",$0}'
done<$ifile \
	| awk '/Clean/'\
	| cut -f 1,2,4,5,10-17 \
	| sed -e 's/(/\t/g' -e 's/)//g' -e 's/%//g' \
	| awk -v OFS="\t" 'BEGIN{print "Sample\tlane\tread pairs\tdata size\teffective(%)\tQ20(%)\tQ30(%)\tGC(%)\terror rate(%)"}{print $1,$2,$3,$4,$5,$6";"$7,$8";"$9,$10";"$11,$12";"$13}'\
	> $ofile
