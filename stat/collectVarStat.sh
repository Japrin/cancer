#!/bin/bash

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <in> <out> [snp/indel]"
	printf "\t<in> format: sample\t.genome.final.vcf.gz.stat.txt file\n"
	exit 1
fi

fList=$1
ofile=$2

varType="snp"
if [ $# -gt 2 ]
then
	varType=$3
fi

oDir=`dirname $ofile`

mkdir -p $oDir
tfile="$oDir/var.stat.Title.col"
printf "sampleID\nexonic\nexonic,splicing\nsplicing\nncRNA_exonic\nncRNA_splicing\nncRNA_UTR3\nncRNA_UTR5\nncRNA_intronic\nUTR5\nUTR3\nUTR5,UTR3\nintronic\nupstream\ndownstream\nupstream,downstream\nintergenic\nTotal\nHet\nHom\nHet_ratio\ntransition\ntransvertion\nts/tv\nin dbSNP (percentage)\nnovel\nnovel ts\nnovel tv\nnovel ts/tv\n" > $tfile
while read ID f
do
	printf "$ID\n" > $oDir/$ID.var.stat
	sed -n  -e '3,3p' $f \
		| transpose.sh \
		| sed '1,1d' \
		>> $oDir/$ID.var.stat
	tfile="$tfile $oDir/$ID.var.stat"
done<$fList
if [ "$varType" == "snp" ]
then
	paste $tfile > $ofile.genome.txt
else
	paste $tfile |  sed -n -e '1,21p' -e '25,26p' > $ofile.genome.txt
fi
rm $tfile

tfile="$oDir/var.stat.Title.col"
printf "sampleID\nframeshift_insertion\nframeshift_deletion\nframeshift_substitution\nstopgain_SNV\nstoploss_SNV\nnonframeshift_insertion\nnonframeshift_deletion\nnonframeshift_substitution\nmissense_SNV\nsynonymous_SNV\nunknown\nTotal\nHet\nHom\nHet_ratio\ntransition\ntransvertion\nts/tv\nin dbSNP (percentage)\nnovel\nnovel ts\nnovel tv\nnovel ts/tv\nNS:S ratio\n" > $tfile
while read ID f
do
	printf "$ID\n" > $oDir/$ID.var.stat
	sed -n  -e '6,6p' $f \
		| transpose.sh \
		| sed '1,1d' \
		>> $oDir/$ID.var.stat
	tfile="$tfile $oDir/$ID.var.stat"
done<$fList
if [ "$varType" == "snp" ]
then
	paste $tfile | sed -e '2,4d' -e '7,9d' > $ofile.exome.txt
else
	paste $tfile | sed -e '10,11d' -e '17,19d' -e '22,24d' > $ofile.exome.txt
fi
rm $tfile

