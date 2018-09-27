#!/bin/bash -eu

echo "*** Recalibrating base quality ***"

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
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
java -Djava.io.tmpdir=$o_dir -Xms${optM}g -Xmx${optM}g -jar $GATK/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-R $REF \
	--knownSites $knownSites \
	--disable_indel_quals \
	-I $f \
	-o ${o/.bam/.grp} \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov ContextCovariate \
	-rf BadCigar
	#-et NO_ET \
	#-K $GATKKey

echo ">> Table recalibration"
if [ "`grep -v '#' ${o/.bam/.grp} | grep -v "EOF" | wc -l`" = "1" ]
then
	cp $f $o
else
	java -Djava.io.tmpdir=$o_dir -Xms${optM}g -Xmx${optM}g -jar $GATK/GenomeAnalysisTK.jar \
	-T PrintReads \
	-R $REF \
	-I $f \
	-o $o \
	-BQSR ${o/.bam/.grp} \
	--disable_indel_quals \
	-EOQ \
	-rf BadCigar
	#-et NO_ET \
	#-K $GATKKey
fi
rm -f ${o/.bam/.bai}
samtools index $o

rm $f

echo "*** Finished recalibrating base quality ***"
