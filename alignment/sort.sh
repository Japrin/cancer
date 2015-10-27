#!/bin/bash


iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
optT=4
optM=16

while getopts c:p:t:m: opt
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
		optT=$OPTARG
		;;
	m)
		optM=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-t threads, default 4] [-m memrory(GB), default 12] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-c iniFile] [-t threads, default 4] [-m memrory(GB), default 12] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi

echo begin at: `date`
echo "*** Aligning reads ***"

source $iniFile

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
lib=$5
lane=$6

mkdir -p $outDir

if [ -f $outDir/$sampleID.$lib.$lane.sort.bam ]
then
    echo "## $outDir/$sampleID.$lib.$lane.bam and $outDir/$sampleID.$lib.$lane.sort.bam exists, skip this step"
	exit 0
fi
## link files ##
sai1="$outDir/`basename $fq1`.sai"

if [ "$fq2" != '-' ]
then
	#ln -f -s $fq2 $outDir/`basename $fq2`
	#fq2=$outDir/`basename $fq2`
	sai2="$outDir/`basename $fq2`.sai"
fi

optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`
echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
java -Xmx${optM}g -jar $picardDIR/SortSam.jar \
		I=$outDir/$sampleID.$lib.$lane.bam \
		O=$outDir/$sampleID.$lib.$lane.sort.bam \
		MAX_RECORDS_IN_RAM=$max_reads \
		TMP_DIR=$outDir \
		SO=coordinate \
		VALIDATION_STRINGENCY=SILENT

#rm $outDir/$sampleID.$lib.$lane.bam

echo end at: `date`
echo "*** Finished aligning reads ***"
