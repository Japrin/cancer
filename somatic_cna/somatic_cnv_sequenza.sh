#!/bin/bash

shDir=`dirname $0`
iniFile="`dirname $0`/../parameter/init_human.hg38.sh"
source $iniFile
_refData=$REF
TR=""
snvfile=""
normalSNPFile=""
gender="F"

while getopts c:b:r:f:g:s: opt
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
		if [ -f $OPTARG ]
		then
			TR=$OPTARG
		else
	        echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] [-f reference] [-g gender,default \"M\"] [-s somatic snv file] <sampleID> <normalBam> <tumorBam> <outDir>"
			exit 1
		fi
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
    s)
		if [ -f $OPTARG ]
		then
			snvfile=$OPTARG
		fi
		;;
	b)
		normalSNPFile=$OPTARG
		;;
	g)
		gender=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] [-f reference] [-g gender,default \"F\"] [-s somatic snv file] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] [-f reference] [-g gender,default \"F\"] [-s somatic snv file] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

##source $iniFile
refData=$_refData
REF=$_refData

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir

### 
if [ ! -f "$outDir/$sampleID.all.binned.seqz.gz" ];then
    
    sequenza-utils bam2seqz -n $normalBam -t $tumorBam --fasta $refData \
        -gc $bundleDir/Homo_sapiens_assembly38.gc50Base.wig.gz \
        --parallel 8 \
        --chromosome "chr"{1..22} chrX \
        -o $outDir/$sampleID.seqz.gz
    
    (
    for xx in {1..22} X
    do
        bgzip -cd $outDir/${sampleID}_chr$xx.seqz.gz
    done
    ) | awk 'NR==1 || !/^chromosome/' | bgzip -c > $outDir/$sampleID.all.seqz.gz
    
    sequenza-utils seqz_binning --seqz $outDir/$sampleID.all.seqz.gz -w 50 -o $outDir/$sampleID.all.binned.seqz.gz 

fi

if [ -f "$snvfile" ] && [ ! -f "$outDir/$sampleID.somatic_snvs.seqz.gz" ];then
    sequenza-utils snp2seqz \
        --vcf $snvfile \
        -gc $bundleDir/Homo_sapiens_assembly38.gc50Base.wig.gz \
        --preset strelka2_som \
        --output $outDir/$sampleID.somatic_snvs.seqz.gz

    sequenza-utils seqz_merge \
        --seqz1 $outDir/$sampleID.somatic_snvs.seqz.gz \
        --seqz2 $outDir/$sampleID.all.binned.seqz.gz \
        --output $outDir/$sampleID.somatic_plus_snps.seqz.gz
fi

if [ ! -f "$outDir/${sampleID}_segments.txt" ];then
    $shDir/somatic_cnv_sequenza.inR.R \
        -i $outDir/$sampleID.all.binned.seqz.gz \
        -o $outDir/R \
        -x $gender \
        -s $sampleID
fi

echo end at: `date`
