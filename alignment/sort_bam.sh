#!/bin/bash


iniFile="`dirname $0`/../parameter/init_human.sh"
optT=4
optM=16
optA="coordinate"

while getopts c:p:t:m:a: opt
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
		optT=$OPTARG
		;;
    a)
        optA=$OPTARG
        ;;
	m)
		optM=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 12] [-a sort by, default \"coordinate\", alt \"queryname\"] inbam <outbam>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 1 ]
then 
	echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 12] [-a sort by, default \"coordinate\", alt \"queryname\"] inbam <outbam>"
	exit 1
fi

echo begin at: `date`
echo "*** Aligning reads ***"

source $iniFile

#mkdir -p $outDir
inbam=$1
if [ $# -lt 2 ];then
    outbam=${inbam/.bam/.sort.bam}
else
    outbam=$2
fi
outDir=`dirname $outbam`

echo begin at: `date`

optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`
echo "... using $max_reads reads in memory (parameter optM: $optM*1.5)"
java -Xmx${optM}g -jar $picardDIR/picard.jar SortSam \
		I=$inbam \
		O=$outbam \
		MAX_RECORDS_IN_RAM=$max_reads \
		TMP_DIR=$outDir \
		SO=$optA \
		VALIDATION_STRINGENCY=SILENT
if [ "$optA" == "coordinate" ];then
    samtools index $outbam
fi
#rm $outDir/$sampleID.$lib.$lane.bam

echo end at: `date`
