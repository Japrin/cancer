#!/bin/bash


iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"
optT=4
optM=12
optI=""
optA=0
optB=1000
_refData="/PROJ/GR/share/medinfo.00database/UCSC/human/chromFaMasked/mrFAST/hg19.masked.fa"
snpfile="/PROJ/GR/share/medinfo.00database/UCSC/human/chromFaMasked/mrsFAST/dbsnp_137.b37.onlySNP.index"
while getopts c:f:s:t:m:i:a:b: opt
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
	f)
		_refData=$OPTARG
		;;
	s)
		snpfile=$OPTARG
		;;
	t)
		optT=$OPTARG
		;;
	m)
		optM=$OPTARG
		;;
	i)
		optI=".$OPTARG"
		;;
	a)
		optA=$OPTARG
		;;
	b)
		optB=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-f refData] [-t threads, default 4] [-m memrory(GB), default 12] [-a min_distance, default 0] [-b max_distance, default 1000] [-i part] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-c iniFile] [-f refData] [-t threads, default 4] [-m memrory(GB), default 12] [-a min_distance, default 0] [-b max_distance, default 1000] [-i part] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi

echo begin at: `date`
echo "*** Aligning reads ***"

source $iniFile
refData=$_refData
REF=$refData

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
#
if [ -f $outDir/$sampleID.$lib.$lane.sort.bam ]
then
    echo "## $outDir/$sampleID.$lib.$lane.sort.bam exists, skip this step"
	exit 0
fi
#
if [ ! -f "$fq1" ]
then
	echo "## $fq1 not exists !!!!"
	exit 0
fi


## bwa ##
if [ ! -f "$outDir/$sampleID.$lib.$lane$optI.bam" ];then
	printf ">>>> runing bwa <<<<\n"
	if [ "$fq2" != '-' ]; then
		mrfast  --search  $refData \
			--pe \
			--seq1 $fq1 \
			--seq2 $fq2 \
			--min $optA \
			--max $optB \
			-o $outDir/$sampleID.$lib.$lane$optI \
			-u $outDir/$sampleID.$lib.$lane$optI.nohit \
			--seqcomp \
			--outcomp \
			--sample $sampleID \
			--rg $sampleID \
			--lib $lib
	else
		mrfast  --search  $refData \
			--seq $fq1 \
			-o $outDir/$sampleID.$lib.$lane$optI\
			-u $outDir/$sampleID.$lib.$lane$optI.nohit \
			--seqcomp \
			--outcomp \
			--sample $sampleID \
			--rg $sampleID \
			--lib $lib
	fi
	gzip -f -v $outDir/$sampleID.$lib.$lane$optI.nohit
	gzip -cd $outDir/$sampleID.$lib.$lane$optI.gz \
		| samtools view -b -S -t $refData.fai - > $outDir/$sampleID.$lib.$lane$optI.bam 
	#java -jar $picardDIR/AddOrReplaceReadGroups.jar \
	#	I=$outDir/$sampleID.$lib.$lane$optI.tmp.bam \
	#	O=$outDir/$sampleID.$lib.$lane$optI.bam \
	#	SO=coordinate \
	#	ID=$sampleID \
	#	LB=$lib \
	#	PL="illumina" \
	#	PU=${lib}_${lane} \
	#	SM=$sampleID \
	#	CN="novogene"
	if [ -f "$outDir/$sampleID.$lib.$lane$optI.bam" ];then
		rm -f $outDir/$sampleID.$lib.$lane$optI.gz
		rm -f $outDir/$sampleID.$lib.$lane$optI.tmp.bam
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
	samtools index $outDir/$sampleID.$lib.$lane$optI.sort.bam
	if [ -f "$outDir/$sampleID.$lib.$lane$optI.sort.bam" ];then
		rm $outDir/$sampleID.$lib.$lane$optI.bam
	fi
fi

echo end at: `date`
echo "*** Finished aligning reads ***"
