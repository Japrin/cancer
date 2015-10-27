#!/bin/bash


iniFile="/PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh"
_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optP=""
optT=4
optM=16
optI=""
optE=50

while getopts c:p:t:m:f:e:i: opt
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
	p)
		optP=$OPTARG
		;;
	t)
		optT=$OPTARG
		;;
	m)
		optM=$OPTARG
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
	e)
		optE=$OPTARG
		;;
	i)
		optI=".$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-p other 'bwa aln' parameters ] [-f reference] [-e gap extension, default 50] [-t threads, default 4] [-m memrory(GB), default 16] [-i part] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-c iniFile] [-p other 'bwa aln' parameters ] [-f reference] [-e gap extension, default 50] [-t threads, default 4] [-m memrory(GB), default 16] [-i part] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi

echo begin at: `date`
echo "*** Aligning reads ***"

source $iniFile
refData=$_refData
REF=$_refData

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
lib=$5
lane=$6

mkdir -p $outDir

if [ -f $outDir/$sampleID.$lib.$lane$optI.sort.bam ]
then
    echo "## $outDir/$sampleID.$lib.$lane$optI.sort.bam exists, skip this step"
	exit 0
fi

if [ -f $outDir/$sampleID.$lib.$lane.sort.bam ]
then
    echo "## $outDir/$sampleID.$lib.$lane.sort.bam exists, skip this step"
	exit 0
fi

if [ ! -f "$fq1" ]
then
	echo "## $fq1 not exists !!!!"
	exit 0
fi
## link files ##
#ln -f -s $fq1 $outDir/`basename $fq1`
#fq1=$outDir/`basename $fq1`
sai1="$outDir/`basename $fq1`.sai"

if [ "$fq2" != '-' ]
then
	#ln -f -s $fq2 $outDir/`basename $fq2`
	#fq2=$outDir/`basename $fq2`
	sai2="$outDir/`basename $fq2`.sai"
fi

## bwa ##
if [ ! -f "$sai1" ];then
	printf ">>>> runing bwa <<<<\n"
	if [ "$fq2" != '-' ]; then
		bwa aln $optP -e $optE -i 15 -q 15 -t $optT $refData $fq1 > $sai1
		bwa aln $optP -e $optE -i 15 -q 15 -t $optT $refData $fq2 > $sai2
	else
		bwa aln $optP -e $optE -i 15 -q 15 -t $optT $refData $fq1 > $sai1
	fi
fi

if [ ! -f "$outDir/$sampleID.$lib.$lane$optI.bam" ];then
	if [ "$fq2" != '-' ]; then
		bwa sampe -P -r "@RG\tID:$sampleID\tSM:$sampleID\tLB:$lib\tPU:${lib}_${lane}\tPL:illumina\tCN:BIOPIC" $refData $sai1 $sai2 $fq1 $fq2 \
			| samtools view -b -S -t $refData.fai - > $outDir/$sampleID.$lib.$lane$optI.bam
	else
		bwa samse -r "@RG\tID:$sampleID\tSM:$sampleID\tLB:$lib\tPU:${lib}_${lane}\tPL:illumina\tCN:BIOPIC" $refData $sai1 $fq1 \
			| samtools view -b -S -t $refData.fai - > $outDir/$sampleID.$lib.$lane$optI.bam
	fi

fi

if [ ! -f "$outDir/$sampleID.$lib.$lane$optI.sort.bam" ];then
	optM=`echo "scale=0;$optM/1.5" | bc`
	max_reads=`echo 250000*$optM | bc`
	echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
	java -Xmx${optM}g -jar $picardDIR/SortSam.jar \
			I=$outDir/$sampleID.$lib.$lane$optI.bam \
			O=$outDir/$sampleID.$lib.$lane$optI.sort.bam \
			MAX_RECORDS_IN_RAM=$max_reads \
			TMP_DIR=$outDir \
			SO=coordinate \
			VALIDATION_STRINGENCY=SILENT
	if [ -f "$outDir/$sampleID.$lib.$lane$optI.sort.bam" ];then
		rm $outDir/$sampleID.$lib.$lane$optI.bam
	fi
fi

echo end at: `date`
echo "*** Finished aligning reads ***"
