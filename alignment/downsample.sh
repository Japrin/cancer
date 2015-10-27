#!/bin/bash -eu

echo "*** downsample bam ***"

iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"
optM=16
optP=0.1

while getopts c:m:p: opt
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
	m)
		optM=$OPTARG
		;;
	p)
		optP=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 6] [-p probability, default 0.1] <in> <out>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 6] [-p probability, default 0.1] <in> <out>"
	exit 1
fi

source $iniFile

#f=`cd \`dirname $1\`; pwd`/`basename $1`
#o=`cd \`dirname $2\`; pwd`/`basename $2`
f=$1
o=$2
o_dir=`dirname $o`

if [ -f $o ]
then
    echo "## $o exists, skip this step"
	exit 0
fi

if [ ! -f $f.bai ]
then
	rm -f ${f/.bam/.bai}
	samtools index $f
fi

echo ">>> Marking duplicates"
optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`
echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
java -Xmx${optM}g -jar $PICARD/DownsampleSam.jar \
	TMP_DIR=$o_dir \
	I=${f} \
	O=${o} \
	P=$optP \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=$max_reads

rm -f ${o/.bam/.bai}
samtools index $o
echo "*** Finished removing duplicates ***"
