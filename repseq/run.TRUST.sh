#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optBam=false
while getopts b opt
do
	case $opt in 
	b)	
		optBam=true
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 <outDir> <sampleID> <inbam>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))


if [ $# -lt 3 ]
then 
	echo "Usage: $0 <outDir> <sampleID> <inbam>"
	exit 1
fi

outDir=$1
sampleID=$2
inbam=$3

optT=4

mkdir -p $outDir

source $iniFile
##. /usr/share/Modules/init/bash
##export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
##module load java/1.8.0_112
##module load bamUtil/1.0.12

echo begin at: `date` `hostname`

#inbam="/WPSnew/zhenglt/work/TCR_chunhong/extra_dataset/TCGA/cart/rnaseq/aln/009a36da-7fed-4c0e-b1ef-d139499fbe42/c3413264-49a1-4968-8acd-43d463b4250f_gdc_realn_rehead.bam"

trustBin="/Share/BP/zhenglt/01.bin/repSeq/trust/SupplementarySoftware/TRUST.pyc"
if [ ! -f "$inbam.bai" ];then
    ln -s ${inbam/.bam/.bai} $inbam.bai
fi
cd $outDir
python $trustBin -f $inbam -a > $sampleID.log 2>&1

echo end at: `date`
