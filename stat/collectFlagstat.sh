#!/bin/bash

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <in> <out>"
	printf "\t<in> format: sample\t.flagstat.txt file\n"
	exit 1
fi

fList=$1
ofile=$2

oDir=`dirname $ofile`
mkdir -p $oDir
tfile="$oDir/Flagstat.Title.col"
printf "Sample\nTotal\nDuplicate\nMapped\nProperly mapped\nPE mapped\nSE mapped\nwith mate mapped to a different chr\nwith mate mapped to a different chr (mapQ>=5)\n" > $tfile
while read ID f
do
	sed -n -e '1,3p' -e '7,11p' $f \
		| awk -v OFS="\t" 'NR==1{t=$1;print "'$ID'"}{printf "%d (%4.2f%)\n",$1,100*$1/t}' \
		> $oDir/$ID.flagstat
	tfile="$tfile $oDir/$ID.flagstat"
done<$fList

paste $tfile > $ofile
rm $tfile
