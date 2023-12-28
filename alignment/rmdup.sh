#!/bin/bash -eu

echo "*** Removing duplicates ***"

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.hg38.sh"
optM=16

while getopts c:m: opt
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
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 6] <in> <out> [remove, default: false]"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 6] <in> <out> [remove, default: false]"
	exit 1
fi

source $iniFile

rmdup="false"
if [ $# -gt 2 ]
then
	rmdup=$3
fi

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=`cd \`dirname $2\`; pwd`/`basename $2`

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

	###M=${o/.bam/.metrics} \
echo ">>> Marking duplicates"
optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`
echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
java -Xmx${optM}g -jar $PICARD MarkDuplicates \
	TMP_DIR=$o_dir \
	I=${f} \
	O=${o} \
	M=${o%.bam}.metrics \
	VALIDATION_STRINGENCY=SILENT \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=$rmdup \
	MAX_RECORDS_IN_RAM=$max_reads

rm -f ${o%.bam}.bai
samtools index $o

rm $f

echo "*** Finished removing duplicates ***"
