#!/bin/bash

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"

optM=10

while getopts c:m: opt
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
	m)
		optM=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 10] <sampleID> <inbam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 10] <sampleID> <inbam> <outDir>"
	exit 1
fi

echo begin at: `date`
source $iniFile
module load RNA-SeQC/1.1.8
#module load bwa/0.6.2
module load bwa/0.7.12

sampleID=$1
inbam=$2
outDir=$3

mkdir -p $outDir

GENE_MODEL="/DBS/DB_temp/zhangLab/ucsc/annotation/hg19/mybuild/hg19_knownGene.gtf"
GCFile="/DBS/DB_temp/zhangLab/ucsc/annotation/hg19/mybuild/hg19_knownGene.Tx.GC"
REF="/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/bwa_0.7.12/hg19.fa"
RRNA_INDEX="/DBS/DB_temp/zhangLab/gencode/release_19/hisat/human_all_rRNA"

if [ -f "$inbam" ];then
	JM=`echo "scale=0;$optM/1.5" | bc`
	java -Xmx${JM}g -jar $RNASeQC \
		-n 1000 \
		-s "$sampleID|$inbam|NA" \
		-t $GENE_MODEL \
		-r $REF \
		-o $outDir \
		-strat gc -gc $GCFile \
		-BWArRNA "$RRNA_INDEX.fa"
fi


echo end at: `date`
