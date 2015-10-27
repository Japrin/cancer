#!/bin/bash

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <in> <out>"
	printf "\t<in> format: sample\t.coverage.bychr.txt file\n"
	exit 1
fi

fList=$1
ofile=$2

ofile1=$ofile.1.txt
ofile2=$ofile.2.txt

rm -f $ofile1 $ofile2
#printf "\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\tX\tY\n" > $ofile1
#printf "\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\tX\tY\n" > $ofile2
while read ID f
do
	cut -f 1,4 $f \
		| transpose.sh \
		| awk -F"\t" -v OFS="\t" 'NR==1{print $0}NR==2{$1="'$ID'";print $0}' >> $ofile1
	cut -f 1,6 $f \
		| transpose.sh \
		| awk -F"\t" -v OFS="\t" 'NR==1{print $0}NR==2{$1="'$ID'";print $0}' >> $ofile2
done<$fList
awk 'NR==1||!/^chr/' $ofile1 > $ofile1.tmp && mv $ofile1.tmp $ofile1
awk 'NR==1||!/^chr/' $ofile2 > $ofile2.tmp && mv $ofile2.tmp $ofile2
plot.coverage.byChr.R $ofile1 $ofile2 $ofile
