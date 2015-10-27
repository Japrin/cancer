#!/bin/bash

echo "*** somatic INDEL by GATK ***"

TR=""
optTR=""
iniFile="/PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh"

while getopts c:r: opt
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
	r)	
		if [ -f $OPTARG ]
		then
			TR="$OPTARG"
			optTR="-L $TR"
		else
			echo "WARNING: invalid target file ($OPTARG), no target will be used"
		fi
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir

java -Xmx4G -jar $gatkJAR \
	-T SomaticIndelDetector \
	-l INFO \
	-R $refData $optTR \
	-I:normal $normalBam \
	-I:tumor $tumorBam \
	-metrics $outDir/$sampleID.SomaticIndelDetector.callability.txt \
	-bed $outDir/$sampleID.SomaticIndelDetector.brief.bed \
	-verbose $outDir/$sampleID.SomaticIndelDetector.detailed.txt \
	--window_size 300 \
	-mnr 10000000 \
	-et NO_ET \
	-K $GATKKey \
	--filter_expressions "T_COV<6||N_COV<4||T_INDEL_F<0.3||T_INDEL_CF<0.7" \
	-o $outDir/$sampleID.SomaticIndelDetector.raw.vcf

GATKSomaticIndel_filter.pl --NDepthTh 6 --TDepthTh 8 --NFreqTh 0.01 --TFreqTh 0.1 --MMTh 2 --NQSMMTh 0.01 $outDir/$sampleID.SomaticIndelDetector.raw.vcf > $outDir/$sampleID.SomaticIndelDetector.call.vcf 2> $outDir/$sampleID.SomaticIndelDetector.call.err
#awk '/^#/|| /SOMATIC/' $outDir/$sampleID.SomaticIndelDetector.raw.vcf > $outDir/$sampleID.SomaticIndelDetector.call.vcf

#ss=`head -1 $outDir/$sampleID.SomaticIndelDetector.call.err`

echo end at: `date`
echo "*** Finished somatic Indel by SomaticIndelDetector ***"
