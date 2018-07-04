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
	    echo "Usage: $0 [-t threads, default 8] <inDir> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 1 ]
then 
	echo "Usage: $0 [-t threads, default 8] <inDir> "
	exit 1
fi

source $iniFile
module load gcc/4.9.2

inDir=$1
TMPDIR="/Share/BP/zhenglt/tmp"
export PATH="/Share/BP/zhenglt/01.bin/RepSeq/tracer":$PATH
TRACER_CONFIG="/Share/BP/zhenglt/01.bin/RepSeq/tracer/tracer.conf"
export IGDATA="/DBS/DB_temp/zhangLab/IMGT/igblast"
#-auxiliary_data "/DBS/DB_temp/zhangLab/IMGT/igblast/optional_file/human_gl.aux"

echo begin at: `date`

echo tracer summarise --keep_inkt -c $TRACER_CONFIG --graph_format "svg" $inDir
tracer summarise --keep_inkt -c $TRACER_CONFIG --graph_format "svg" $inDir

echo end at: `date`
