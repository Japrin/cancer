#!/bin/bash

echo "***  ***"

iniFile="/WPS/GR/zhengliangtao/02pipeline/novo.med.seq/cancer/parameter/init_human.sh"

while getopts c: opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid reference file ($OPTARG), default will be used"
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] <sampleID> <inFile> <outFile>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] <sampleID> <inFile> <outFile>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
infile=$2
outfile=$3

mkdir -p `dirname $outfile`
cd `dirname $outfile`

R CMD BATCH --save --restore "--args infile=\"$infile\" outfile=\"$outfile\"" /WPS/GR/zhengliangtao/01bin/novoHumanReseq_v1.0.1/app/novoHumanReseq/bin/cancer/scna_CNAnorm.R

echo end at: `date`
echo "*** Finished exome CNV ***"
