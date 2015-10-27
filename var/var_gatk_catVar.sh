#!/bin/bash -eu
echo "*** variant calling by GATK ***"

iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"

optM=6

while getopts c:r:p:m: opt
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
	r)	
		TR="$OPTARG"
		optTR="-L $TR"
		;;
	m)
		optM=$OPTARG
		;;
	p)
		opt_P="$OPTARG";
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 6] <sampleID> <byChrDir> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-m memrory(GB), default 6] <sampleID> <byChrDir> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile 

sampleID=$1
byChrDir=$2
outDir=$3
vcfList="$byChrDir/vcf.raw.list"
ofile="$outDir/$sampleID.GATK.var.raw.vcf"

for i in {1..22} X Y MT
do 
	printf "$byChrDir/$sampleID.$i.GATK.var.raw.vcf\n"
done  > $vcfList


optV=`awk '{print "-V "$1}' $vcfList | xargs`
optOther=""
optOther="-assumeSorted"

if [ ! -f "$ofile" ];then
	java -Xmx${optM}g -Djava.io.tmpdir=$outDir -cp $GATK/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.CatVariants \
		-R $REF \
	 	$optV \
		-out $ofile $optOther
else
	echo "$ofile exists, skip this step"
fi

echo end at: `date`
