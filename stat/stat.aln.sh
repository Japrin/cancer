#!/bin/bash

echo "*** bam stat ***"

TR=""
optTR=""
iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/bwa_0.7.12/human_g1k_v37_decoy.fasta"

while getopts c:r:f: opt
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
			TR="$OPTARG"
			optTR="-l $TR"
		else
			echo "WARNING: invalid target file ($OPTARG), no target will be used"
		fi
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-f reference] <bam> <outDir> <ID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-f reference] <bam> <outDir> <ID>"
	exit 1
fi

echo begin at: `date`

source $iniFile 

module load samtools/0.1.19

refData=$_refData
REF=$_refData

bam=$1
outDir=$2
sampleID=$3

if [ ! -f $bam.bai ]
then
	rm -f ${bam/.bam/.bai}
	samtools index $bam
fi

mkdir -p $outDir
if [ -f "$TR" ]
then
	echo ">> Target Region coverage stat"
	stat_exom.pl -i $bam -r $TR -o $outDir -plot
	cal.coverage.sh -c $iniFile -r $TR -f $REF $sampleID $bam $outDir
else
	echo ">> Whole genome coverage stat"
	depth.addn.pl -q 0 -Q 0 -n $PIPELINE/cancer/stat/b37.NBlock.larger20bp.bed $bam $outDir > $outDir/$sampleID.cov.txt
	cal.coverage.sh -c $iniFile -r $PIPELINE/cancer/stat/b37.chr25Region.bed -f $REF $sampleID $bam $outDir
fi
samtools flagstat $bam > $outDir/$sampleID.stat.aln.flagstat.txt

##
tfile="$outDir/$sampleID.Flagstat.Title.col"
printf "Sample\nTotal\nDuplicate\nMapped\nProperly mapped\nPE mapped\nSE mapped\nwith mate mapped to a different chr\nwith mate mapped to a different chr (mapQ>=5)\n" > $tfile
sed -n -e '1,3p' -e '7,11p' $outDir/$sampleID.stat.aln.flagstat.txt | awk -v OFS="\t" 'NR==1{t=$1;print "'$sampleID'"}{printf "%d (%4.2f%)\n",$1,100*$1/t}' > $outDir/tmp123456
tfile="$tfile $outDir/tmp123456"
paste $tfile > $outDir/tmp1234567
rm $tfile
tfile=$outDir/tmp1234567
if [ -f "$TR" ]
then
	cat $tfile $outDir/information.xlsx > $outDir/${sampleID}_mapping_coverage.txt
else
	cat $tfile $outDir/$sampleID.cov.txt > $outDir/${sampleID}_mapping_coverage.txt
fi
rm $tfile

echo end at: `date`
echo "*** Finished bam stat ***"
