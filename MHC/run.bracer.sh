#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optT=8
optS="Hsap"

while getopts t:s: opt
do
	case $opt in 
	t)
		optT=$OPTARG
		;;
    s)
        optS=$OPTARG
        ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
        echo "Usage: $0 [-t threads, default 8] [-s Hsap (default) or Mmus] <outDir> <fq1> <fq2> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
    echo "Usage: $0 [-t threads, default 8] [-s Hsap (default) or Mmus] <outDir> <fq1> <fq2> <sampleID>"
	exit 1
fi

source $iniFile
module load gcc/4.9.2
module unload python/2.7.8
module load python/3.6.1
module unload java/1.7.0_79
module load java/1.8.0_144

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
TMPDIR="/Share/BP/zhenglt/tmp"
#export PATH="/Share/BP/zhenglt/01.bin/RepSeq/tracer":$PATH
BRACER_CONFIG="/Share/BP/zhenglt/01.bin/RepSeq/bracer/bracer.conf"
export IGDATA="/DBS/DB_temp/zhangLab/IMGT/igblast"
#-auxiliary_data "/DBS/DB_temp/zhangLab/IMGT/igblast/optional_file/human_gl.aux"

####export LC_CTYPE="en_US.UTF-8" 

### test 
bracer test -p $optT -c $BRACER_CONFIG --infer_lineage

echo begin at: `date`

#mkdir -p $outDir

echo end at: `date`
