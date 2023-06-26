#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"
#_refData="/WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
###_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optT="--localcores 4"
optM="--localmem 16"
optS="human"
optA="TCR"
optP=""

while getopts c:t:m:s:a:p: opt
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
		optT="--localcores $OPTARG"
		;;
	m)
		optM="--localmem $OPTARG"
		;;
    s)
        optS=$OPTARG
        ;;
    a)
        optA=$OPTARG
        ;;
    p)
        optP=$OPTARG
        ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
        echo "Usage: $0 [-c iniFile] [-t threads, default 4] [-m memrory(GB), default 16] [-s species, default human] [-a vdj type, TCR(default) or BCR] [-p cellranger parameters] <outDir> <inDir> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
    echo "Usage: $0 [-c iniFile] [-t threads, default 4] [-m memrory(GB), default 16] [-s species, default human] [-a vdj type, TCR(default) or BCR] [-p cellranger parameters] <outDir> <inDir> <sampleID>"
	exit 1
fi

echo begin at: `date`

#source $iniFile

outDir=$1
inDir=$2
sampleID=$3

if [ "$optS" == "human" ];then
    #transcriptomeDir="/WPSnew/zhenglt/00.database/ensemble/10X/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0"
    transcriptomeDir="/workspace/zhengliangtao/00.database/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
elif [ "$optS" == "mouse" ];then
    transcriptomeDir="/WPSnew/zhenglt/00.database/ensemble/10X/refdata-cellranger-vdj-GRCm38-alts-ensembl-2.2.0"
fi

mkdir -p $outDir
###export MODULESHOME=/usr/share/Modules
###. /usr/share/Modules/init/bash
###export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
#module load cellranger/2.1.1
#module load cellranger/3.0.0
module load cellranger/5.0.0

echo begin at: `date`
cd $outDir

newInDirVDJ=$outDir/$optA.FQ.$sampleID
mkdir -p $newInDirVDJ
ln -s `ls $inDir/*_R1*.gz` $newInDirVDJ/${sampleID}_S1_L001_R1_001.fastq.gz
ln -s `ls $inDir/*_R2*.gz` $newInDirVDJ/${sampleID}_S1_L001_R2_001.fastq.gz

#ln -s $inDir/*_R1*.gz $newInDirVDJ/
#ln -s $inDir/*_R2*.gz $newInDirVDJ/

echo cellranger vdj --id=$sampleID $optT $optM \
    --fastqs=$newInDirVDJ \
    --sample=$sampleID \
    --reference=$transcriptomeDir $optP

cellranger vdj --id=$sampleID $optT $optM \
    --fastqs=$newInDirVDJ \
    --sample=$sampleID \
    --reference=$transcriptomeDir $optP

    ####--expect-cells=3000
echo end at: `date`



