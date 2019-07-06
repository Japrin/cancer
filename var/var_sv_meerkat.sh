#!/bin/bash -eu


iniFile="`dirname $0`/../parameter/init_human.sh"

optN=""

while getopts c:n: opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid file ($OPTARG), default will be used"
		fi
		;;
	n)	
		optN="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG"
		echo "Usage: $0 [-c iniFile] [-n normal result] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-n normal result] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile
#REF="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/hg19/hg19.fa"
meerkatDir="/WPSnew/zhenglt/01.bin/cancer/sv/meerkat/Meerkat.0.189"

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

## redefin path
export LD_LIBRARY_PATH=$meerkatDir/lib:$LD_LIBRARY_PATH
export PATH=/WPSnew/zhenglt/01.bin/var/samtools/samtools-0.1.19/mybuild/bin:/WPSnew/zhenglt/01.bin/aln/bwa/bwa-0.6.2:/opt/bio/ncbi/bin:/WPSnew/zhenglt/01.bin/cancer/blat/blatSrc35/mybuild/x86_64:$PATH

mkdir -p $outDir/normal $outDir/tumor
#ln -s $normalBam $outDir/normal/$sampleID.normal.bam
#ln -s $normalBam.bai $outDir/normal/$sampleID.normal.bam.bai
#ln -s $tumorBam $outDir/tumor/$sampleID.tumor.bam
#ln -s $tumorBam.bai $outDir/tumor/$sampleID.tumor.bam.bai

##perl $meerkatDir/scripts/pre_process.pl -b $outDir/tumor/$sampleID.tumor.bam -I $REF -A $REF.fai -t 8 -k 2000 -c 45 -s 20
##perl $meerkatDir/scripts/meerkat.pl     -b $outDir/tumor/$sampleID.tumor.bam -F `dirname $REF`/BioDB -t 8 -d 5 -s 20 -Q 10 -p 5 -o 3
###perl $meerkatDir/scripts/mechanism.pl   -b $outDir/tumor/$sampleID.tumor.bam -R $HumanDB/rmsk.txt

#perl $meerkatDir/scripts/pre_process.pl -b $outDir/tumor/$sampleID.tumor.bam -I $REF -A $REF.fai -t 8 -k 2000 -c 45 -s 20 -l 0 && \
perl $meerkatDir/scripts/meerkat.pl -b $outDir/tumor/$sampleID.tumor.bam -F `dirname $REF`/BioDB -t 8 -d 5 -s 20 -Q 10 -p 4 -o 2 -l 0 -m 0 && \
perl $meerkatDir/scripts/mechanism.pl   -b $outDir/tumor/$sampleID.tumor.bam -R $HumanDB/hg19_rmsk.txt.gz &

#if [ ! -f "$optN" ]
#then
#	##perl $meerkatDir/scripts/pre_process.pl -b $outDir/normal/$sampleID.normal.bam -I $REF -A $REF.fai -t 8 -k 2000 -c 45 -s 20
#	##perl $meerkatDir/scripts/meerkat.pl     -b $outDir/normal/$sampleID.normal.bam -F `dirname $REF`/BioDB -t 8 -d 5 -s 20 -Q 10 -p 5 -o 3
#	###perl $meerkatDir/scripts/mechanism.pl   -b $outDir/normal/$sampleID.normal.bam -R $HumanDB/rmsk.txt
#	perl $meerkatDir/scripts/pre_process.pl -b $outDir/normal/$sampleID.normal.bam -I $REF -A $REF.fai -t 8 -k 2000 -c 45 -s 20 -l 0 && \
#	perl $meerkatDir/scripts/meerkat.pl -b $outDir/normal/$sampleID.normal.bam -F `dirname $REF`/BioDB -t 8 -d 5 -s 20 -Q 10 -p 4 -o 2 -l 0 -m 0 && \
#	perl $meerkatDir/scripts/mechanism.pl   -b $outDir/normal/$sampleID.normal.bam -R $HumanDB/hg19_rmsk.txt.gz &
#fi

wait


### filter

somatica=$outDir/tumor/$sampleID.tumor.somatica.variants
somaticb=$outDir/tumor/$sampleID.tumor.somaticb.variants
somaticc=$outDir/tumor/$sampleID.tumor.somaticc.variants
somaticd=$outDir/tumor/$sampleID.tumor.somaticd.variants
somatice=$outDir/tumor/$sampleID.tumor.somatice.variants
somaticf=$outDir/tumor/$sampleID.tumor.somaticf.variants
somaticg=$outDir/tumor/$sampleID.tumor.somaticg.variants
normal_blacklistrg=""
tumor_blacklistrg=""
normal_bam=$outDir/normal/$sampleID.normal.bam
tumor_bam=$outDir/tumor/$sampleID.tumor.bam

perl $meerkatDir/scripts/somatic_sv.pl -i $outDir/tumor/$sampleID.tumor.variants -o $somatica -F $outDir/normal -l 1000 -R $HumanDB/hg19_rmsk.txt.gz
perl $meerkatDir/scripts/somatic_sv.pl -i $somatica -o $somaticb -n 1 -Q 10 -B $normal_bam -I $outDir/normal/$sampleID.normal.isinfo -D 5 -K $outDir/normal/$sampleID.normal.blacklist.gz  -R $HumanDB/hg19_rmsk.txt.gz
perl $meerkatDir/scripts/somatic_sv.pl -i $somaticb -o $somaticc -u 1 -Q 10 -B $normal_bam  -R $HumanDB/hg19_rmsk.txt.gz -K $outDir/normal/$sampleID.normal.blacklist.gz
perl $meerkatDir/scripts/somatic_sv.pl -i $somaticc -o $somaticd -f 1 -Q 10 -B $normal_bam  -R $HumanDB/hg19_rmsk.txt.gz -K $outDir/normal/$sampleID.normal.blacklist.gz
perl $meerkatDir/scripts/somatic_sv.pl -i $somaticd -o $somatice -e 1 -Q 10 -B $tumor_bam -I $outDir/tumor/$sampleID.tumor.isinfo -D 5 -R $HumanDB/hg19_rmsk.txt.gz -K $outDir/normal/$sampleID.normal.blacklist.gz
perl $meerkatDir/scripts/somatic_sv.pl -i $somatice -o $somaticf -z 1 -R $HumanDB/hg19_rmsk.txt.gz
perl $meerkatDir/scripts/somatic_sv.pl -i $somaticf -o $somaticg -d 40 -t 20 -R $HumanDB/hg19_rmsk.txt.gz

## gff format
perl $PIPELINE/cancer/var/var_sv_meerkat2gff.pl -p $sampleID $somaticg > $somaticg.gff
$PIPELINE/cancer/all/var_annotation.sh $somaticg.gff $sampleID

## vcf format
outputVCF=$outDir/tumor/$sampleID.tumor.somatic.vcf
perl $meerkatDir/scripts/meerkat2vcf.pl -i $somaticg -o $outputVCF -F `dirname $REF`/BioDB -H $meerkatDir/dataset/headerfile

## fusion
perl $meerkatDir/scripts/fusions.pl -i $somaticg -G $HumanDB/hg19_refGene_sorted.txt

echo end at: `date`
