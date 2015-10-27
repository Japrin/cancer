#!/bin/bash

shDir=`dirname $0`
iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
TR=""

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
			TR=$OPTARG
		else
			echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] <sampleID> <normalBam> <tumorBam> <outDir>"
			exit 1
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir

#samtools index $normalBam
#samtools index $tumorBam

if [ ! -f "$outDir/$sampleID.normal.sample_interval_summary" ];then
	java -Djava.io.tmpdir=$outDir -Xmx5g -jar $gatkJAR -T DepthOfCoverage \
		-mmq 1 \
		-mbq 13 \
		-omitLocusTable \
		-omitBaseOutput \
		-R $REF \
		-I $normalBam \
		-L $TR \
		-o $outDir/$sampleID.normal \
		-et NO_ET \
		-K $gatkKey
	sed -i 's/^chr//' $outDir/$sampleID.normal.sample_interval_summary
fi

if [ ! -f "$outDir/$sampleID.tumor.sample_interval_summary" ];then
	java -Djava.io.tmpdir=$outDir -Xmx5g -jar $gatkJAR -T DepthOfCoverage \
		-mmq 1 \
		-mbq 13 \
		-omitLocusTable \
		-omitBaseOutput \
		-R $REF \
		-I $tumorBam \
		-L $TR \
		-o $outDir/$sampleID.tumor \
		-et NO_ET \
		-K $gatkKey
	sed -i 's/^chr//' $outDir/$sampleID.tumor.sample_interval_summary
fi

#Rscript $shDir/ecnv_analyze_cnv.R $outDir/$sampleID.normal.sample_interval_summary $outDir/$sampleID.tumor.sample_interval_summary $outDir/$sampleID

(
awk 'NR==1{print $0"\ttarget.number"}' $outDir/$sampleID.cnv.cnv.txt
sed '1,1d' $outDir/$sampleID.cnv.cnv.txt | intersectBed -a stdin -b $outDir/$sampleID.cnv.exon.lrr.txt -c
) > $outDir/$sampleID.cnv.cnv.targetN.txt

awk -F"\t" 'NR==1 || ($7!=2 && $13>=10)' $outDir/$sampleID.cnv.cnv.targetN.txt > $outDir/$sampleID.cnv.cnv.targetN.filter.txt

Rscript $shDir/plot_ExomeCNV.R $sampleID $outDir/$sampleID.cnv.exon.lrr.txt $outDir/$sampleID.cnv.cnv.targetN.txt $outDir

echo end at: `date`
