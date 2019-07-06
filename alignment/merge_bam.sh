#!/bin/bash

echo "*** Merge BAM ***"

iniFile="`dirname $0`/../parameter/init_human.sh"
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
        echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 16] <output> <bam>..."
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then
        echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 16] <output> <bam>..."
        exit 1
fi

source $iniFile

mkdir -p `dirname $1`
out=`cd \`dirname $1\`; pwd`/`basename $1`
shift 1

if [ -f $out ]
then
    echo "## $out exists, skip this step"
	exit 0
fi

bams=''
bams2=''

o_dir=`dirname $out`

for f in $*
do
	if [ -f $f ];then
		#f=`cd \`dirname $f\`; pwd`/`basename $f`
		bams="$bams I=$f"
		bams2="$bams2 $f"
	fi
done

echo ">>> Merging BAMs into $out"
if [ $# -gt 1 ]
then
	###samtools merge $out $bams
	optM=`echo "scale=0;$optM/1.5" | bc`
	max_reads=`echo 250000*$optM | bc`
	echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
	java -Xmx${optM}g -jar $picardDIR/picard.jar MergeSamFiles \
		$bams \
		O=$out \
		TMP_DIR=$o_dir \
		SO=coordinate \
		MSD=true \
		MAX_RECORDS_IN_RAM=$max_reads \
		VALIDATION_STRINGENCY=SILENT
	#rm $bams2
else
	#cp $bams2 $out
	ln -f $bams2 $out
fi

samtools index $out
echo "*** Finished Merge BAM ***"
