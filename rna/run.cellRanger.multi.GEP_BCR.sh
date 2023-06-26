#!/bin/bash

shDir=`dirname $0`
#iniFile="$sDir/../parameter/init_human.sh"
templateFile="$shDir/template.cellRanger.multi.config.GEP_BCR.csv"
###_refData="/WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
###_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optT="16"
optM="128"
optS="human"
optP=""
optMode=""

while getopts c:a:t:m:s:p: opt
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
	a)
		optMode="$OPTARG"
		;;
	t)
		optT="$OPTARG"
		;;
	m)
		optM="$OPTARG"
		;;
    s)
        optS=$OPTARG
        ;;
    p)
        optP=$OPTARG
        ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-a mode. one of 'fastp', and ''. default ''] [-t threads, default 16] [-m memrory(GB), default 128] [-p cellranger parameters] <outDir> <sampleID> <inDir_GEP> <inDir_BCR>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-a mode. one of 'fastp', and ''. default ''] [-t threads, default 16] [-m memrory(GB), default 128] [-p cellranger parameters] <outDir> <sampleID> <inDir_GEP> <inDir_BCR>"
	exit 1
fi

echo begin at: `date`

#source $iniFile

outDir=$1
sampleID=$2
inDirGEP=$3
inDirBCR=$4

mkdir -p $outDir
module load cellranger/5.0.0
module load fastp/0.23.2

#outDir=`pwd`/OUT.cellranger
#mkdir -p $outDir
#mkdir -p $outDir/fq

echo begin at: `date`
cd $outDir

rm -r $sampleID

newInDirGEP=$outDir/GEP.FQ.$sampleID
newInDirBCR=$outDir/BCR.FQ.$sampleID
mkdir -p $newInDirGEP
mkdir -p $newInDirBCR

if [[ $optMode == "fastp" ]];then
    ### [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
    fastp -i `ls $inDirGEP/*_R1*.gz` -I `ls $inDirGEP/*_R2*.gz` \
          -o $newInDirGEP/${sampleID}-GEP_S1_L001_R1_001.fastq.gz -O $newInDirGEP/${sampleID}-GEP_S1_L001_R2_001.fastq.gz \
          --length_required 75 --detect_adapter_for_pe \
          --thread $optT
    fastp -i `ls $inDirBCR/*_R1*.gz` -I `ls $inDirBCR/*_R2*.gz` \
          -o $newInDirBCR/${sampleID}-BCR_S1_L001_R1_001.fastq.gz -O $newInDirBCR/${sampleID}-BCR_S1_L001_R2_001.fastq.gz \
          --length_required 75 --detect_adapter_for_pe \
          --thread $optT
    echo fastp done.
else
#    ln -s $inDirGEP/*_R1*.gz $newInDirGEP/
#    ln -s $inDirGEP/*_R2*.gz $newInDirGEP/
#    ln -s $inDirBCR/*_R1*.gz $newInDirBCR/
#    ln -s $inDirBCR/*_R2*.gz $newInDirBCR/
    ln -s `ls $inDirGEP/*_R1*.gz` $newInDirGEP/${sampleID}-GEP_S1_L001_R1_001.fastq.gz
    ln -s `ls $inDirGEP/*_R2*.gz` $newInDirGEP/${sampleID}-GEP_S1_L001_R2_001.fastq.gz
    ln -s `ls $inDirBCR/*_R1*.gz` $newInDirBCR/${sampleID}-BCR_S1_L001_R1_001.fastq.gz
    ln -s `ls $inDirBCR/*_R2*.gz` $newInDirBCR/${sampleID}-BCR_S1_L001_R2_001.fastq.gz
fi
inDirGEP=$newInDirGEP
inDirBCR=$newInDirBCR

sed -e "s#SampleID_GEP#${sampleID}-GEP#" \
    -e "s#SampleID_BCR#${sampleID}-BCR#" \
    -e "s#FQ_DIR_GEP#${inDirGEP}#" \
    -e "s#FQ_DIR_BCR#${inDirBCR}#" \
    $templateFile \
    > $outDir/cellRanger.multi.config.$sampleID.csv

cellranger multi \
    --id=$sampleID \
    --csv=$outDir/cellRanger.multi.config.$sampleID.csv \
    --localmem=$optM --localcores=$optT $optP

echo end at: `date`



