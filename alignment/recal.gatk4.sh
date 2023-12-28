#!/bin/bash -eu

echo "*** Recalibrating base quality ***"

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.hg38.sh"
optM=8
optNT=1
optNCT=4

while getopts c: opt
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
		echo "Usage: $0 [-c iniFile] <in> <out>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile] <in> <out>"
	exit 1
fi

source $iniFile

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

echo ">>> Counting covariates"
optM=`echo "scale=0;$optM/1.2" | bc`
gatk --java-options "-Djava.io.tmpdir=$o_dir -Xms${optM}g -Xmx${optM}g -Dsamjdk.compression_level=5" BaseRecalibrator \
	-R $REF \
	--known-sites $knownSites \
	-I $f \
	-O ${o/.bam/.grp}

echo ">> Table recalibration"
if [ "`grep -v '#' ${o/.bam/.grp} | grep -v "EOF" | wc -l`" = "1" ]
then
	ln $f $o
else
	gatk --java-options "-Djava.io.tmpdir=$o_dir -Xms${optM}g -Xmx${optM}g -Dsamjdk.compression_level=5" ApplyBQSR \
	    -R $REF \
        -I $f \
        -O $o \
        --emit-original-quals \
        --bqsr-recal-file ${o/.bam/.grp}
fi
rm -f ${o/.bam/.bai}
samtools index $o

rm $f $f.bai

echo "*** Finished recalibrating base quality ***"
