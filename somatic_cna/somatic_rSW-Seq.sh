#!/bin/bash

echo "*** somatic SNV by varScan ***"

NCount=1
TCount=1
SWScore=80
iniFile="/WPS/GR/zhengliangtao/02pipeline/novo.med.seq/cancer/parameter/init_human.sh"

while getopts c:N:T:S: opt
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
	N)	
		NCount=$OPTARG
		;;
	T)
		TCount=$OPTARG
		;;
	S)
		SWScore=$OPTARG
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-N NCount(1)] [-T TCount(1)] [-S SW score(80)] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-N NCount(1)] [-T TCount(1)] [-S SW score(80)] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

for i in {1..22} X Y
do
	printf "process chr$i\n";
	mkdir -p $outDir/chr$i

	normalTxt=$outDir/chr$i/$sampleID.rSW-seq.normal.input
	tumorTxt=$outDir/chr$i/$sampleID.rSW-seq.tumor.input
	normalNFile=$outDir/chr$i/$sampleID.rSW-seq.normal.count
	tumorNFile=$outDir/chr$i/$sampleID.rSW-seq.tumor.count

	samtools view $normalBam -F 0x404 -q 37 chr$i  | awk -F"\t" '{sum++;print $4 > "'$normalTxt'"} END{print sum > "'$normalNFile'"}'
	samtools view $tumorBam -F 0x404 -q 37 chr$i  | awk -F"\t" '{sum++;print $4 > "'$tumorTxt'"} END{print sum > "'$tumorNFile'"}'
done
NCount=`cat $outDir/chr*/$sampleID.rSW-seq.normal.count | awk '{sum+=$1}END{print sum}'`
TCount=`cat $outDir/chr*/$sampleID.rSW-seq.tumor.count | awk '{sum+=$1}END{print sum}'`
printf "" > $outDir/$sampleID.rSW-Seq.bed
for i in {1..22} X Y
do
	normalTxt=$outDir/chr$i/$sampleID.rSW-seq.normal.input
	tumorTxt=$outDir/chr$i/$sampleID.rSW-seq.tumor.input
	/WPS/GR/zhengliangtao/01bin/rSW-seq/rSW-seq \
		$tumorTxt \
		$normalTxt \
		$TCount \
		$NCount \
		$SWScore \
	> $outDir/chr$i/$sampleID.rSW-seq.txt
	awk -v OFS="\t" '/^gain/||/^loss/ {if($1~/loss/){ t=$4;$4=$5;$5=t} r=$4/$5; print "chr'$i'",$2,$3,$3-$2,$1,r,$4,$5,$6,$7;}' $outDir/chr$i/$sampleID.rSW-seq.txt \
		| sort -k 1,1 -k 2g,2 -k 3g,3 \
		>> $outDir/$sampleID.rSW-Seq.bed
done

awk '$9>=80 && ($6>1.4 || $6<0.6)' $outDir/$sampleID.rSW-Seq.bed > $outDir/$sampleID.rSW-Seq.SWScore80.bed

windowBed -a $outDir/$sampleID.rSW-Seq.SWScore80.bed \
	-b /WPS/GR/zhengliangtao/00database/human/GATK_bundle/1.5/hg19/addn.all.more1000.bed \
	-w 50000 \
	-v \
	> $outDir/$sampleID.rSW-Seq.final.bed

echo end at: `date`
echo "*** Finished somatic SNV by varScan ***"
