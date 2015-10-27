#!/bin/bash


iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"

binSize=30000

while getopts b:c: opt
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
	b)
		binSize="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-b binSize, default 30000] <ID> <outDir> <bamlist file>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-b binSize, default 30000] <ID> <outDir> <bamlist file>"
	exit 1
fi

echo begin at: `date`

source $iniFile 

ID=$1
outDir=$2
bamlistfile=$3

mkdir -p $outDir

echo begin at: `date`
cd $outDir

Rscript $PIPELINE/var/run.cn.mops.R -b $binSize -i $bamlistfile -o $outDir/$ID


echo end at: `date`
