#!/bin/bash

echo "***  ***"

iniFile="/WPS/GR/zhengliangtao/02pipeline/novo.med.seq/cancer/parameter/init_human.sh"
optP=""
optW=""

while getopts c:p:w: opt
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
	p)
		optP="--p-value $OPTARG"
		;;
	w)
		optW="--window-size $OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-p p-value(default 0.001)] [-w window-size(default determined by log2 and p-value)] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-p p-value(default 0.001)] [-w window-size(default determined by log2 and p-value)] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir

cd $outDir

samtools view -F 0x404 -q 37 $normalBam | perl -lane 'if($F[2]=~/^(chr)*[0-9XY]+$/){ print "$F[2]\t$F[3]"}' > $outDir/$sampleID.normal.hits
samtools view -F 0x404 -q 37 $tumorBam  | perl -lane 'if($F[2]=~/^(chr)*[0-9XY]+$/){ print "$F[2]\t$F[3]"}' > $outDir/$sampleID.tumor.hits

/WPS/GR/zhengliangtao/01bin/CNV-seq/cnv-seq/cnv-seq.pl \
	--test $outDir/$sampleID.tumor.hits \
	--ref  $outDir/$sampleID.normal.hits \
	--genome human \
	--Rexe /usr/local/bin/R $optP $optW

head -1 $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv \
	> $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted
sed '1,1d' $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv > $outDir/.tmp
vcf_sort.pl $outDir/.tmp >> $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted
rm $outDir/.tmp

/WPS/GR/zhengliangtao/01bin/novoHumanReseq_v1.0.1/app/novoHumanReseq/bin/cancer/CNVseq.print.CNV.R \
	$outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted \
	| sed 's/chrchr/chr/' \
	> $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment 
head -1 $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment > $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment.sorted
sed -e '1,1d' $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment | sort -k 2,2 -k 3g,3 -k 4g,4 > $outDir/.tmp
vcf_sort.pl -i 1 $outDir/.tmp >> $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment.sorted
rm $outDir/.tmp

#other format
awk -v OFS="\t" 'NR>1 {f=($6>0)?"gain":"loss";print $2,$3,$4,$5,f,$6,"NA","NA",$7,$1}' $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment.sorted > $outDir/$sampleID.CNV-Seq.bed

windowBed -a $outDir/$sampleID.CNV-Seq.bed \
	-b /WPS/GR/zhengliangtao/00database/human/GATK_bundle/1.5/hg19/addn.all.more1000.bed \
	-w 50000 \
	-v \
	> $outDir/$sampleID.CNV-Seq.final.bed 

## plot (paper)
CNVseq.plot.paper.R $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted
## plot (refGene track)
#mkdir -p $outDir/CNVR
#sed 's/CNVR_//' $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment.sorted | awk -v OFS="\t" 'NR==1 {print "chr","id"} NR>1{print $2,$1}' > $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment.sorted.id.list
#CNVseq.plot.R $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted $outDir/$sampleID.tumor.hits-vs-$sampleID.normal.hits.log2-0.6.pvalue-0.001.minw-4.cnv.sorted.segment.sorted.id.list $outDir/CNVR

echo end at: `date`
echo "*** Finished ***"
