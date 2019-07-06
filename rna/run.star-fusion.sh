#!/bin/bash

#!/bin/bash -eu

iniFile="`dirname $0`/../parameter/init_human.sh"
optNT=8

while getopts c:p: opt
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
	p)
		optNT=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-p numOfCPU, default 8] <outDir> <fq1> <fq2> <aid>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-p numOfCPU, default 8] <outDir> <fq1> <fq2> <aid>"
	exit 1
fi

source $iniFile

STAR_FUSION_DB_Dir=/WPSnew/zhenglt/00.database/tools/star-fusion/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir

module load blast/2.8.1+
module load STAR/2.6.1d
module load STAR-Fusion/1.5.0

outDir=$1
fq1=$2
fq2=$3
aid=$4

mkdir -p $outDir

echo begin at: `date`

STAR-Fusion --genome_lib_dir $STAR_FUSION_DB_Dir \
             --left_fq $fq1 \
             --right_fq $fq2 \
             --CPU $optNT \
             --output_dir $outDir

echo end at: `date`

