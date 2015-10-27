#!/bin/bash


shDir=`dirname $0`
iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
gender="M"

while getopts c:g: opt
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
		gender=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-g gender, default \"M\" for male, alter is \"F\" for female ] <snp (vcf.gz)> <normalBam> <tumorBam> <outDir> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 5 ]
then 
	echo "Usage: $0 [-c iniFile] [-g gender, default \"M\" for male, alter is \"F\" for female ] <snp (vcf.gz)> <normalBam> <tumorBam> <outDir> <sampleID>"
	exit 1
fi

echo begin at: `date`

source $iniFile

snpFile=$1
normalBam=$2
tumorBam=$3
outDir=$4
sampleID=$5

mkdir -p $outDir

if [ ! -f "$outDir/$sampleID.hetSNV.RC" ]; then
	bgzip -cd $snpFile \
		| /Share/BP/zhenglt/02.pipeline/cancer/rna/mbased.filterSNV.pl -n -f -a 10 -b 0.1 -g $gender \
		> $outDir/$sampleID.hetSNV.bed
	samtools mpileup -f $refData -l $outDir/$sampleID.hetSNV.bed $normalBam $tumorBam > $outDir/$sampleID.pileup
	countMultiPileup.pl -maxNormal 1 -minDepth 1 $outDir/$sampleID.pileup > $outDir/$sampleID.hetSNV.RC
fi

awk -F"\t" -v OFS="\t" 'BEGIN{print "chr\tposition\tcoverage\tbaf"} {print "chr"$1,$2,$4+$5,$5 }' $outDir/$sampleID.hetSNV.RC > $outDir/$sampleID.normal.LOH.input
awk -F"\t" -v OFS="\t" 'BEGIN{print "chr\tposition\tcoverage\tbaf"} {print "chr"$1,$2,$6+$7,$7 }' $outDir/$sampleID.hetSNV.RC > $outDir/$sampleID.tumor.LOH.input

Rscript $shDir/ecnv_analyze_loh.R $outDir/$sampleID.normal.LOH.input $outDir/$sampleID.tumor.LOH.input $outDir/$sampleID

awk -F"\t" -v OFS="\t" '{$3=$3+1; print $0}' $outDir/$sampleID.the.loh.txt \
	| intersectBed -a stdin -b $outDir/$sampleID.eLOH.loh.txt -wo \
	| $shDir/ExomeCNV.LOH.postProcess.pl \
	> $outDir/$sampleID.ExomeCNV.LOH

awk 'NR==1 || ($8=="TRUE" && $10 > 5)' $outDir/$sampleID.ExomeCNV.LOH > $outDir/$sampleID.ExomeCNV.LOH.txt

echo end at: `date`
