#!/bin/bash -eu

echo "*** Realigning targeted regions ***"

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
optM=8
optNT=1
optNCT=4
optL=''

while getopts c:m:l: opt
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
    l)
		if [ -f $OPTARG ]
		then
			optL="-L $OPTARG"
        fi
        ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-l interval] <sampleID> <inbam> <outDir> <REF>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-l interval] <sampleID> <inbam> <outDir> <REF>"
	exit 1
fi

source $iniFile

module unload java/1.7.0_79
module load java/1.8.0_144
module unload gatk/3.3-0
module load gatk/3.8-0
### with ERCC
#export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/withERCC.v1/hg19.order.after.gsnap.fa
### without ERCC
#export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/hg19.order.after.gsnap.fa
knownSites=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/dbsnp_138.hg19.reorder.vcf

sampleID=$1
inBam=$2
outDir=$3
REF=$4

mkdir -p $outDir


optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`

java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $REF \
    -I $inBam $optL \
    -dontUseSoftClippedBases \
    -stand_call_conf 20.0 \
	-nt $optNT \
    --num_cpu_threads_per_data_thread 4 \
    -o $outDir/$sampleID.GATK.RNA.recal.var.vcf

java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $REF \
    -V $outDir/$sampleID.GATK.RNA.recal.var.vcf \
    -window 35 -cluster 3 \
    -filterName FS -filter "FS > 30.0" \
    -filterName QD -filter "QD < 2.0" \
    -o $outDir/$sampleID.GATK.RNA.recal.var.flt.vcf


echo end at: `date` 
