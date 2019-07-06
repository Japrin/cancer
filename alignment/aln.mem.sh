#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"
###iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
#_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/bwa_0.7.12//human_g1k_v37_decoy.fasta"
optT=4
optM=16
optI=""
optP=""

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
	p)
		optP=$OPTARG
		;;
	i)
		optI=".$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-f reference] [-p extra parameters (Phred64) ] [-t threads, default 4] [-m memrory(GB), default 16] [-i part] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-c iniFile] [-f reference] [-p extra parameters (Phred64) ] [-t threads, default 4] [-m memrory(GB), default 16] [-i part] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi

echo begin at: `date`
echo "*** Aligning reads ***"

source $iniFile
#refData=$_refData
#REF=$_refData

module load seqtk/1.3

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
lib=$5
lane=$6
BinDir=`dirname $0`

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

#ifq1=$outDir/$sampleID.input.R1.fq.gz
#ifq2=$outDir/$sampleID.input.R2.fq.gz
ifq1=$outDir/$sampleID.$lib.$lane$optI.R1.fq.gz
ifq2=$outDir/$sampleID.$lib.$lane$optI.R2.fq.gz

if [ ! -f "$outDir/$sampleID.$lib.$lane$optI.bam" ];then
	BQE=`$BinDir/check.BQ.64.33.pl $fq1`
	echo "detect base quality encoding is in Phred$BQE format !"
	
	if [ "$fq2" != '-' ]; then
		if [ "$BQE" == "64" ];then
			echo "convert input fq encoding in Phred64 format to encoding in Phred33 format ......"
			mkfifo $ifq1
			mkfifo $ifq2
			seqtk seq -VQ64 $fq1  > $ifq1 &
			seqtk seq -VQ64 $fq2  > $ifq2 &
			#seqtk seq -VQ64 $fq1 | gzip -c -f > $ifq1
			#seqtk seq -VQ64 $fq2 | gzip -c -f > $ifq2
		else
			echo "make link file: $ifq1 and $ifq2"
			ln -s -f $fq1 $ifq1
			ln -s -f $fq2 $ifq2
		fi

		bwa mem -M -t $optT -R "@RG\tID:$sampleID\tSM:$sampleID\tLB:$lib\tPU:${lane}\tPL:illumina" $refData $ifq1 $ifq2 \
			| samtools view -bS - > $outDir/$sampleID.$lib.$lane$optI.bam &
		wait
		#rm $ifq1 $ifq2

	else
		if [ "$BQE" == "64" ];then
			echo "convert input fq encoding in Phred64 format to encoding in Phred33 format ......"
			seqtk seq -VQ64 $fq1 | gzip -c -f > $ifq1
		else
			ln -s $fq1 $ifq1
		fi
		bwa mem -M -t $optT -R "@RG\tID:$sampleID\tSM:$sampleID\tLB:$lib\tPU:${lane}PL:illumina" $refData $ifq1 \
			| samtools view -bS - > $outDir/$sampleID.$lib.$lane$optI.bam
		rm $ifq1
	fi

fi

if [ ! -f "$outDir/$sampleID.$lib.$lane$optI.sort.bam" ];then
	optM=`echo "scale=0;$optM/1.5" | bc`
	max_reads=`echo 250000*$optM | bc`
	echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
	java -Xmx${optM}g -jar $picardDIR/picard.jar SortSam \
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
