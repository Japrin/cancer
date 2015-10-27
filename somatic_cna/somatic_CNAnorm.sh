#!/bin/bash

echo "***  ***"

winSize=1000
gcFile="/PROJ/GR/share/medinfo.00database/CNAnorm/gc1000Base.txt.gz"
iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"

while getopts c:w:g: opt
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
	w)	
		winSize=$OPTARG
		;;
	g)
		if [ -f $OPTARG ]
		then
			gcFile="$OPTARG"
		else
			echo "WARNING: invalid gc file ($OPTARG), default will be used"
		fi
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-w winSize(1000)] [-g gcFile] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-w winSize(1000)] [-g gcFile] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir

/PROJ/GR/share/Software/medinfo/01bin/CNAnorm/bam2windows.pl \
	$tumorBam \
	$normalBam \
	-V \
	-w $winSize \
	-gc $gcFile \
	-q 37 \
	-d $outDir \
	-st $sampleID.CNAnorm.tumor.temp  \
	-sc $sampleID.CNAnorm.normal.temp \
	-ts -cs \
	> $outDir/$sampleID.CNAnorm.input

R CMD BATCH --save --restore "--args infile=\"$outDir/$sampleID.CNAnorm.input\" outfile=\"$outDir/$sampleID.CNAnorm.output\"" $PIPELINE/somatic_cna/scna_CNAnorm.R
$PIPELINE/somatic_cna/CNAnorm2Bed.pl -w $winSize $outDir/$sampleID.CNAnorm.output | awk '$9!="NA" && ($9<1.5||$9>2.5)' > $outDir/$sampleID.CNAnorm.bed
windowBed -a $outDir/$sampleID.CNAnorm.bed -b /PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/all.NBlock.larger1000bp.bed -w 50000 -v > $outDir/$sampleID.CNAnorm.final.bed

echo end at: `date`
echo "*** Finished ***"
