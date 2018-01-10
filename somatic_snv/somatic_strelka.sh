#!/bin/bash

echo "*** somatic mutation by strelka ***"

TR=""
optTR=""
iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
#_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optZ="N"
optG=""

shDir=`dirname $0`

while getopts c:r:f:z:g: opt
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
	g)
		if [ -f $OPTARG ];then
			optG=" -g $OPTARG "
		fi
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
	z)
		optZ="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-g gender file] [-f reference] [-z N|Y] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-g gender file] [-f reference] [-z N|Y] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile
#refData=$_refData
#REF=$_refData

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

strelkaDir="/Share/BP/zhenglt/01.bin/strelka/strelka_workflow-1.0.14/mybuild"

#mkdir -p $outDir
config_file=""

if [ -f "$TR" ]; then
	config_file="$strelkaDir/etc/strelka_config_bwa_default_TR.ini"
else
	config_file="$strelkaDir/etc/strelka_config_bwa_default_WGS.ini"
fi
### call ###
if [ ! -f "$outDir/results/passed.somatic.snvs.vcf" ]; then

	echo "good"
	$strelkaDir/bin/configureStrelkaWorkflow.pl \
		--tumor=$tumorBam \
		--normal=$normalBam \
		--ref=$refData \
		--config=$config_file \
		--output-dir=$outDir
	
	make -C $outDir -j 8

fi
if [ -f "$TR" ]; then
	(
		awk '/^#/' $outDir/results/passed.somatic.snvs.vcf
		intersectBed -a $outDir/results/passed.somatic.snvs.vcf   -b $TR | uniq 
	) | $shDir/strelka.addFreq.pl > $outDir/results/$sampleID.strelka.somatic.snvs.filter.vcf

	(
		awk '/^#/' $outDir/results/passed.somatic.indels.vcf
		intersectBed -a $outDir/results/passed.somatic.indels.vcf -b $TR | uniq 
	) | $shDir/strelka.addFreq.pl > $outDir/results/$sampleID.strelka.somatic.indels.filter.vcf
else
	awk '/^#/ || $1~/^([0-9]+|X|Y)$/' $outDir/results/passed.somatic.snvs.vcf   | $shDir/strelka.addFreq.pl > $outDir/results/$sampleID.strelka.somatic.snvs.filter.vcf
	awk '/^#/ || $1~/^([0-9]+|X|Y)$/' $outDir/results/passed.somatic.indels.vcf | $shDir/strelka.addFreq.pl > $outDir/results/$sampleID.strelka.somatic.indels.filter.vcf
fi

if [ "$optZ" = "Y" ] || [ "$optZ" = "y" ]
then
	var_annotation.sh -c $iniFile -a somatic -b SNP   -m strelka.snv   $optG $outDir/results/$sampleID.strelka.somatic.snvs.filter.vcf    N,T
	var_annotation.sh -c $iniFile -a somatic -b INDEL -m strelka.indel $optG $outDir/results/$sampleID.strelka.somatic.indels.filter.vcf  N,T

	### to maf format
	bgzip -cd $outDir/results/$sampleID.strelka.somatic.snvs.filter.reformated.vcf.gz \
		| $PIPELINE/cancer/var/vcf2mafv1.31.strelka.py -t - -s "WES" -a $sampleID \
			-g $HumanDB/hg19_knownGene.SymbolToLocus.txt \
			-o $outDir/results/$sampleID.strelka.somatic.snvs.filter.reformated.maf

	bgzip -cd $outDir/results/$sampleID.strelka.somatic.indels.filter.reformated.vcf.gz \
		| $PIPELINE/cancer/var/vcf2mafv1.31.strelka.py -t - -s "WES" -a $sampleID \
			-g $HumanDB/hg19_knownGene.SymbolToLocus.txt \
			-o $outDir/results/$sampleID.strelka.somatic.indels.filter.reformated.maf
fi

echo end at: `date`
echo "*** Finished somatic mutation by strelka ***"
