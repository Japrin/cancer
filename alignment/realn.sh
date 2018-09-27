#!/bin/bash -eu

echo "*** Realigning targeted regions ***"

iniFile="`dirname $0`/../parameter/init_human.sh"
optM=8
optNT=1
optNCT=4

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

optL=''

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


#### if BQ encoding is Phred+64, use  -fixMisencodedQuals
echo ">>> Determining (small) suspicious intervals which are likely in need of realignment"
optM=`echo "scale=0;$optM/1.2" | bc`
java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$o_dir -jar $GATK/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-I $f \
	-R $REF \
	-o ${o/.bam/.intervals} $optL \
	-nt $optNT \
	-rf BadCigar
	#-et NO_ET \
	#-K $GATKKey

echo ">>> Running the realigner over the targeted intervals"
java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$o_dir -jar $GATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-I $f \
	-R $REF \
	-o $o \
	-targetIntervals ${o/.bam/.intervals} \
	-LOD 5 $optL \
	-nt $optNT \
	-rf BadCigar
	#-et NO_ET \
	#-K $GATKKey

rm -f ${o/.bam/.bai}
samtools index $o

rm $f
#### not nessassary ####
####echo ">>> Fixing the mate pairs and order of the realigned reads"
####java -Xms5g -Xmx5g -jar $PICARD/FixMateInformation.jar \
####	TMP_DIR=$o_dir \
####	INPUT=$o \
####	VALIDATION_STRINGENCY=SILENT \
####	SORT_ORDER=coordinate
####rm -f ${o/.bam/.bai}
####samtools index $o


echo "*** Finished realigning targeted regions ***"
