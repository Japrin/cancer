#!/bin/bash

echo "*** downsample bam ***"

iniFile="/WPS/GR/zhengliangtao/01bin/novoHumanReseq_v1.0.1/app/novoHumanReseq/bin/cancer/init_human.sh"
optA=""

while getopts c:a: opt
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
	a)	
		optA="-a $OPTARG"
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile] [-a downsample to 1/a] <inbam> <outbam>"
	exit 1
fi

inbam=$1
outbam=$2

echo begin at: `date`

source $iniFile

samtools view -F 0x404 -q 30 $inbam | downsample.bam.pl -f | samtools view -Sbt $refData.fai - > $outbam
samtools index $outbam

echo end at: `date`
