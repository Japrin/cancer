#!/bin/bash


iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"

opt_t=8
opt_B=4
opt_n=8
opt_m=1

while getopts c:t:B:n:m: opt
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
            opt_t=$OPTARG
            ;;
        B)
            opt_B=$OPTARG
            ;;
        n)
            opt_n=$OPTARG
            ;;
        m)
            opt_m=$OPTARG
            ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-m mode(1)] [-t threads (8)] [-B gsnap mode (4)] [-n num cores (8)] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-m mode(1)] [-t threads (8)] [-B gsnap mode (4)] [-n num cores (8)] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi

echo begin at: `date`

source $iniFile

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
lib=$5
lane=$6

mkdir -p $outDir

binDir=`dirname $0`

#

echo $binDir/HTSeqGenie.core.R $outDir $sampleID $fq1 $fq2 -m $opt_m -t $opt_t -B $opt_B -n $opt_n -a
$binDir/HTSeqGenie.core.R $outDir $sampleID $fq1 $fq2 -m $opt_m -t $opt_t -B $opt_B -n $opt_n -a


echo end at: `date`

