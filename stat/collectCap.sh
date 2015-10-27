#!/bin/bash

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <in> <out>"
	printf "\t<in> format: sample\tinformation.xlsx file\n"
	exit 1
fi

fList=$1
ofile=$2

oDir=`dirname $ofile`

tfile="$oDir/cap.Title.col"

printf "Sample\nInitial bases on target\nInitial bases near target\nInitial bases on or near target\nTotal effective reads\nTotal effective yield(Mb)\nAverage read length(bp)\nEffective sequences on target(Mb)\nEffective sequences near target(Mb)\nEffective sequences on or near target(Mb)\nFraction of effective bases on target\nFraction of effective bases on or near target\nAverage sequencing depth on target\nAverage sequencing depth near target\nMismatch rate in target region\nMismatch rate in all effective sequence\nBase covered on target\nCoverage of target region\nBase covered near target\nCoverage of flanking region\nFraction of target covered with at least 20x\nFraction of target covered with at least 10x\nFraction of target covered with at least 4x\nFraction of flanking region covered with at least 20x\nFraction of flanking region covered with at least 10x\nFraction of flanking region covered with at least 4x\n" > $tfile

while read ID f
do
	cut -f 2 $f | awk 'BEGIN{print "'$ID'"}{print $0}' > $oDir/$ID.cap
	tfile="$tfile $oDir/$ID.cap"
done < $fList
paste $tfile > $ofile
rm $tfile
