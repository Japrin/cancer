#!/bin/bash

while getopts c: opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid configure file ($OPTARG), default will be used"
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] <pepFile> <output.prefix> <proj> <sampleID> <hla> <len>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-c iniFile] <pepFile> <output.prefix> <proj> <sampleID> <hla> <len>"
	exit 1
fi

pepFile=$1
outPrefix=$2
proj=$3
sampleID=$4
hla=$5
len=$6

#pepFile="/WPS1/zhenglt/work/MHC/ICGC/step0.prepare/OUT/BRCA-US/BRCA-US.ccd4a24b-d8cc-4686-9dee-c98b0c5a8d21/BRCA-US.ccd4a24b-d8cc-4686-9dee-c98b0c5a8d21.annovar.all.k9.exonic_variant_function.codingSeq.peptide.fa"
#outPrefix="BRCA-US.ccd4a24b-d8cc-4686-9dee-c98b0c5a8d21"
#proj="BRCA"
#sampleID="ccd4a24b-d8cc-4686-9dee-c98b0c5a8d21"
#hla="HLA-A*02:05"
#len=9

outDir=`dirname $outPrefix`
mkdir -p $outDir
echo begin at: `date`

/WPS1/zhenglt/work/MHC/ICGC/bin/fa2paddingPeptideFa.pl \
	-o $outPrefix \
	-s $proj.$sampleID \
	$pepFile

source /WPS1/zhenglt/work/MHC/ICGC/bin/MHC.pip.func.sh
runIEDB $sampleID $outPrefix.fa $outDir $hla $len

/WPS1/zhenglt/work/MHC/ICGC/bin/formatIEDBResult.pl \
	-b $outPrefix.bed \
	-p $proj  \
	-s $sampleID \
	$outDir/IEDB.$sampleID.$hla.$len.out \
	> $outDir/IEDB.$sampleID.$hla.$len.out.txt

echo end at: `date`
