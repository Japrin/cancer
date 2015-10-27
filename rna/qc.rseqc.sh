#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optT=8

while getopts c:t: opt
do
	case $opt in 
	t)
		optT=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-t threads, default 8] <outDir> <inbam> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-t threads, default 8] <outDir> <inbam> <sampleID>"
	exit 1
fi

source $iniFile

CDIR=`dirname $0`

outDir=$1
inbam=$2
sampleID=$3

REF_GENE_MODEL_BED="/DBS/DB_temp/zhangLab/ucsc/annotation/hg19/mybuild/hg19_knownGene.bed"
REF_FA="/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/hg19.fa"
CHROM_SIZE_TXT="/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/hg19.chr.length"

mkdir -p $outDir
echo begin at: `date`

geneBody_coverage.py   -r $REF_GENE_MODEL_BED -i $inbam  -o $outDir/$sampleID

inner_distance.py      -r $REF_GENE_MODEL_BED -i $inbam -o $outDir/$sampleID

RPKM_saturation.py     -r $REF_GENE_MODEL_BED -i $inbam -o $outDir/$sampleID -c 1

junction_annotation.py -r $REF_GENE_MODEL_BED -i $inbam -o $outDir/$sampleID

junction_saturation.py -r $REF_GENE_MODEL_BED -i $inbam -o $outDir/$sampleID

read_GC.py                                    -i $inbam -o $outDir/$sampleID

tin.py                 -r $REF_GENE_MODEL_BED -i $inbam > $outDir/$sampleID.tin.out

RNA_fragment_size.py   -r $REF_GENE_MODEL_BED -i $inbam > $outDir/$sampleID.fragSize

$CDIR/geneNumber_saturation.R -i $inbam -o $outDir/$sampleID.geneNumber.saturation -v -s $sampleID \
                --gfeature /DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/withERCC/hg19.knownGene.RData \
                -a 0.02 -b 0.02

gzip -vf $outDir/$sampleID.*.txt
gzip -vf $outDir/$sampleID.*.xls
mv $sampleID.analyzed.tin.xls $outDir/
mv $sampleID.analyzed.summary.txt $outDir/

echo end at: `date`
