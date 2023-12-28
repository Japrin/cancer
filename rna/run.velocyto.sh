#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.hg38.sh"
source $iniFile
_refData=$refData
optT="4"
optM="10"
#optS="human"
optP=""

while getopts c:t:m: opt
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
		optT="$OPTARG"
		;;
	m)
		optM="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-t threads, default 4] [-m memrory(GB), default 10] <sampleID> <outDir> <bamFile> <BCFile>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-t threads, default 4] [-m memrory(GB), default 10] <sampleID> <outDir> <bamFile> <BCFile>"
	exit 1
fi

echo begin at: `date`

sampleID=$1
outDir=$2
bamFile=$3
BCFile=$4

mkdir -p $outDir

##bc="${dir}/filtered_feature_bc_matrix/barcodes.tsv.gz"
samtools --version

#samtools sort -m ${optM}G -t CB -O BAM -@ $optT -o $outDir/cellsorted_${sampleID}.bam $bamFile
#ln -s $bamFile $outDir/${sampleID}.bam

gunzip -c $BCFile > $outDir/${sampleID}.bc

velocyto run \
    -b $outDir/${sampleID}.bc \
    -o $outDir \
    -m /workspace/zhengliangtao/00.database/cellranger/refdata-gex-GRCh38-2020-A/genes/hg38_rmsk.gtf \
    -@ $optT \
    -v \
    $bamFile \
    /workspace/zhengliangtao/00.database/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf

    ##$outDir/${sampleID}.bam 
echo end at: `date`
