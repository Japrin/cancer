#!/bin/bash

echo "***  ***"

iniFile="/WPS/GR/zhengliangtao/02pipeline/novo.med.seq/cancer/parameter/init_human.sh"
bedOPT=""
TRFile=""

while getopts c:b:r: opt
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
	b)
		bedOPT="--bed"
		;;
	r)
		if [ -f $OPTARG ]
		then
			TRFile="$OPTARG"
		else
			echo "WARNING: invalid TR file ($OPTARG)"
			echo "Usage: $0 [-c iniFile] [-b contral is bed file] [-r target region (required)] <sampleID> <control> <test> <outDir>"
			exit 1
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-b contral is bed file] [-r target region (required)] <sampleID> <control> <test> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-b contral is bed file] [-r target region (required)] <sampleID> <control> <test> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p `dirname $outDir` #(contra.py will do this)

/WPS/GR/zhengliangtao/01bin/CONTRA/CONTRA.V2.0/CONTRA.v2.0.3/contra.py \
	--target $TRFile \
	--test $tumorBam \
	--control $normalBam \
	-f $refData \
	-o $outDir \
	--sampleName $sampleID \
	--nomultimapped \
	--plot \
	--largeDeletion $bedOPT

plot_CONTRA.R $outDir/table/$sampleID.CNATable.10rd.10bases.20bins.txt $outDir/table/$sampleID.CNATable.10rd.10bases.20bins.txt.LargeDeletion.txt $outDir/plot

echo end at: `date`
echo "*** Finished ***"
