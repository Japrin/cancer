#!/bin/bash

echo "*** somatic SNV by muTect ***"

TR=""
optTR=""
iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
_refData="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optG=""

while getopts c:r:f:g: opt
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
			optTR="-l $TR"
		else
			TR="$OPTARG"
			optTR="-r $TR"
		fi
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
	g)
		if [ -f $OPTARG ];then
			optG=" -g $OPTARG "
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-g gender file] [-f reference] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-g gender file] [-f reference] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile
refData=$_refData
REF=$_refData

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

muTectDir="/Share/BP/zhenglt/01.bin/mutect/current"

mkdir -p $outDir


### call ###
if [ ! -f "$outDir/$sampleID.muTect.call.vcf" ]; then
	if [ -f "$TR" ]; then
		awk '{print $1":"$2"-"$3}' $TR > $outDir/TR.list
		optTR="--intervals $outDir/TR.list"
	else
		optTR=""
	fi
	/Share/BP/zhenglt/01.bin/java/jre1.6.0_45/bin/java -Xmx6G -Xmx6G -Djava.io.tmpdir=$outDir -jar $muTectDir/muTect.jar \
		-T MuTect \
		-rf BadCigar \
		-dt NONE \
		--reference_sequence $refData \
		--cosmic $muTectDir/b37_cosmic_v54_120711.vcf \
		--dbsnp $muTectDir/dbsnp_132_b37.leftAligned.vcf $optTR \
		--input_file:normal $normalBam \
		--input_file:tumor $tumorBam \
		--vcf $outDir/$sampleID.muTect.call.vcf \
		--out $outDir/$sampleID.muTect.call_stats.out \
		--coverage_file $outDir/$sampleID.coverage.wig.txt
fi
### filter
awk '/^#/||$7!~/REJECT/' $outDir/$sampleID.muTect.call.vcf \
	| awk '/^#/||$1~/^(chr)?([0-9]+|X|Y)/' \
	| $PIPELINE/cancer/var/vcf.SampleReorder.pl \
	> $outDir/$sampleID.muTect.call.filter.vcf

var_annotation.sh -c $iniFile -a somatic -b SNP -m mutect.snv $optG $outDir/$sampleID.muTect.call.filter.vcf  N,T

echo end at: `date`
echo "*** Finished somatic SNV by muTect ***"
