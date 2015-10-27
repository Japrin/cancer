#!/bin/bash


if [ $# -lt 5 ]
then 
	echo "Usage: $0 <infile:*.cnvnator.RD.gz> <chr> <beg> <end> <prefix> <flank>"
	exit 1
fi

infile=$1
chr=$2
beg=$3
end=$4
prefix=$5

flank=50000

if [ $# -gt 5 ]
then
	flank=$6
fi

ll=`echo "$beg-$flank" | bc`
rr=`echo "$end+$flank" | bc`

tmpDir=.
tmpFile=$prefix.for.cnvnator.plot
/PROJ/GR/share/medinfo.02pipeline/bin/tabix $infile "$chr:$ll-$rr" -s 1 -b 2 -e 3 > $tmpFile
/PROJ/GR/share/medinfo.02pipeline/cancer/var/plot.cnvnator.R $tmpFile $prefix $beg $end

awk -F"\t" -v OFS="\t" 'BEGIN{printf "browser full altGraph\ntrack type=bedGraph name= description= visibility=full color=77,110,0 altColor=95,152,205 priority=20 windowingFunction=maximum smoothingWindow=2 alwaysZero=on autoScale=on";} {print $1,$2,$3,log($5)/log(10)}' $tmpFile > $tmpFile.bed


#rm $tmpFile
