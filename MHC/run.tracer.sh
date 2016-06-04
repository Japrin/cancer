#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optT=8

while getopts t: opt
do
	case $opt in 
	t)
		optT=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-t threads, default 8] <outDir> <fq1> <fq2> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-t threads, default 8] <outDir> <fq1> <fq2> <sampleID>"
	exit 1
fi

source $iniFile
module load gcc/4.9.2

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
TMPDIR="/Share/BP/zhenglt/tmp"
export PATH="/Share/BP/zhenglt/01.bin/RepSeq/tracer":$PATH
TRACER_CONFIG="/Share/BP/zhenglt/01.bin/RepSeq/tracer/tracer.conf"
export IGDATA="/DBS/DB_temp/zhangLab/IMGT/igblast"
#-auxiliary_data "/DBS/DB_temp/zhangLab/IMGT/igblast/optional_file/human_gl.aux"

echo begin at: `date`

mkdir -p $outDir
tracer assemble -c $TRACER_CONFIG -p $optT -s Hsap $fq1 $fq2 $sampleID $outDir
#tracer assemble --resume_with_existing_files -c $TRACER_CONFIG -p $optT -s Hsap $fq1 $fq2 $sampleID $outDir

echo end at: `date`
