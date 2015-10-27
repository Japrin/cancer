#!/bin/bash -eu

echo "*** Calling SV using Read-Pair Mapping: ***"

iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"

optG="M"
optP="1"

while getopts c:g:p: opt
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
	g)
		optG="$OPTARG"
		;;
	p)
		optP="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG"
		echo "Usage: $0 [-c iniFile] [-g gender default M] [-p pair type:1 for paired-end 0 for mate-pair, default 1] <sampleID> <outDir> <bam1(tumor)> [bam2(normal) ...]"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-g gender default M] [-p pair type:1 for paired-end 0 for mate-pair, default 1] <sampleID> <outDir> <bam1(tumor)> [bam2(normal) ...]"
	exit 1
fi

echo begin at: `date`

source $iniFile
SV_DETECT_ROOT=/PROJ/GR/share/Software/medinfo/01bin/svdetect/SVDetect_r0.8b

sampleID=$1
outDir=$2
tumorBam=$3
if [ $# -gt 3 ]
then
	normalBam=$4
fi

mkdir -p $outDir

cd $outDir

#perl  $SV_DETECT_ROOT/scripts/BAM_preprocessingPairs.pl -t 1 -p $optP -n 1000000 -s 0 -S 10000 -f 3 -d -o $outDir $tumorBam > $outDir/$sampleID.tumor.preprocessing.log 2>$outDir/$sampleID.tumor.preprocessing.err
#if [ $# -gt 3 ]
#then
#	perl  $SV_DETECT_ROOT/scripts/BAM_preprocessingPairs.pl -t 1 -p $optP -n 1000000 -s 0 -S 10000 -f 3 -d -o $outDir $normalBam > $outDir/$sampleID.normal.preprocessing.log 2>$outDir/$sampleID.normal.preprocessing.err
#fi

#Generation and filtering of links from the sample data
$SV_DETECT_ROOT/bin/SVDetect linking filtering -conf $outDir/sample.sv.conf

if [ $# -gt 3 ]
then
	#Generation and filtering of links from the reference data
	$SV_DETECT_ROOT/bin/SVDetect linking filtering -conf $outDir/reference.sv.conf
	#Comparison of links between the two datasets
	$SV_DETECT_ROOT/bin/SVDetect links2compare -conf $outDir/sample.sv.conf
fi

$SV_DETECT_ROOT/bin/SVDetect links2SV -conf $outDir/sample.sv.conf
###Calculation of depth-of-coverage log-ratios
###$SV_DETECT_ROOT/bin/SVDetect cnv ratio2circos ratio2bedgraph -conf sample.cnv.conf

##Visualization of filtered links and copy-number profiles in Circos
#/PROJ/GR/share/Software/medinfo/01bin/circos/circos-0.64/bin/circos -conf $outDir/sample.circos.conf

echo end at: `date`
