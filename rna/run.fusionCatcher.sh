#!/bin/bash

opt_t=12

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"

while getopts c:t: opt
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
	t)
		opt_t=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-t threads, default 12] <outDir> <sampleID> <tumorFq, \",\"seperated, eg, a.R1.fq,a.R2.fq,b.R1,fq,b.R2.fq> <normalFq or \"-\">"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-t threads, default 12] <outDir> <sampleID> <tumorFq, \",\"seperated, eg, a.R1.fq,a.R2.fq,b.R1,fq,b.R2.fq> <normalFq or \"-\">"
	exit 1
fi

source $iniFile

outDir=$1
sampleID=$2
tumorFq=$3
normalFq=$4

echo begin at: `date`

mkdir -p $outDir

fusionCatcherDir="/Share/BP/zhenglt/01.bin/fusioncatcher/fusioncatcher"

if [ "$normalFq" != '-' ]
then
	$fusionCatcherDir/bin/fusioncatcher \
		--data=$fusionCatcherDir/data/current \
		--input=$tumorFq \
		--normal=$normalFq \
		--output=$outDir \
		--visualization-sam --skip-update-check \
		--threads=$opt_t
else
	$fusionCatcherDir/bin/fusioncatcher \
		--data=$fusionCatcherDir/data/current \
		--input=$tumorFq \
		--output=$outDir \
		--visualization-sam --skip-update-check \
		--threads=$opt_t

fi
echo end at: `date`
