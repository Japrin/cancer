#!/bin/bash -eu

echo "*** Splitting BAM by chromosome ***"

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
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
		echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 16] <chr> <output> <bam>..."
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then
		echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 16] <chr> <output> <bam>..."
        exit 1
fi

source $iniFile 

chr=$1
#out=`cd \`dirname $2\`; pwd`/`basename $2`
out=$2
shift 2

bams=''
bams2=''

o_dir=`dirname $out`
mkdir -p $o_dir

if [ -f $out ]
then
	printf "file exist: $out\n"
    exit 0
fi

for f in $*
do
	f=`cd \`dirname $f\`; pwd`/`basename $f`

	echo ">>> Extracting $chr from BAM: $f"
	
	f_base=`basename $f`
	o=$o_dir/${f_base/.bam/}.$chr.bam
	#o=${f/.bam/}.$chr.bam

	if [ $chr = 'UNK' -o $chr = 'chrU' -o $chr = 'U' ]
	then
		samtools view $f -f 12 -bo $o
	else
		samtools view $f $chr -bo $o
	fi
	bams="$bams I=$o"
	bams2="$bams2 $o"
done

if [ -f $out ]
then
	printf "file exist: $out\n"
    exit 0
fi


echo ">>> Merging $chr BAMs into $out"
if [ $# -gt 1 ]
then
	#samtools merge $out $bams
	optM=`echo "scale=0;$optM/1.2" | bc`
	max_reads=`echo 250000*$optM | bc`
	java -Xmx${optM}g -jar $PICARD/picard.jar MergeSamFiles \
		$bams \
		O=$out \
		TMP_DIR=$o_dir \
		SO=coordinate \
		MSD=true \
		MAX_RECORDS_IN_RAM=$max_reads \
		VALIDATION_STRINGENCY=SILENT
	if [ -f $out ];then
		rm $bams2
	fi
else
	mv $bams2 $out
fi

echo "*** Finished splitting BAM by chromosome ***"
