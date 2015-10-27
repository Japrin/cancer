#!/bin/bash

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
TR=""
optTR=""

while getopts c:r: opt
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
			optTR="-l $TR"
		else
			echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] <sampleID> <normalBam> <tumorBam> <outDir>"
			exit 1
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

shDir=`dirname $0`

mkdir -p $outDir

#if [ ! -f "$outDir/$sampleID.dataRatio.txt" ];then
#	N=`samtools view -q 20 -F 0x400 $normalBam | wc -l`
#	T=`samtools view -q 20 -F 0x400 $tumorBam  | wc -l`
#	dataRatio=`perl -e 'print '$N'/'$T';'`
#	printf "dataRatio\t$dataRatio\n" > $outDir/$sampleID.dataRatio.txt
#else
#	dataRatio=`cut -f 2 $outDir/$sampleID.dataRatio.txt`
#fi
#
#samtoolsDir="/Share/BP/zhenglt/01.bin/samtools/samtools-0.1.18/mybuild/bin"
#
#if [ ! -f "$outDir/$sampleID.somaticCNA.varscan.copynumber.called" ];then
#	printf "No exists: $outDir/$sampleID.somaticCNA.varscan.copynumber.called\n"
#	$samtoolsDir/samtools mpileup -q 1 $optTR -f $refData $normalBam $tumorBam \
#		| perl -F"\t" -ane 'chomp @F;if($F[0]=~/^(\d+|[XY])$/){print}' \
#	    | java -Xmx5G -jar $varScanDIR/VarScan.jar copynumber --mpileup $outDir/$sampleID.somaticCNA.varscan \
#			--min-base-qual 20 \
#			--min-map-qual 20 \
#			--min-coverage 10 \
#			--min-segment-size 10 \
#			--max-segment-size 100 \
#			--p-value 0.01 \
#			--data-ratio $dataRatio
#	
#	java -Xmx5G -jar $varScanDIR/VarScan.jar copyCaller $outDir/$sampleID.somaticCNA.varscan.copynumber \
#		--output-file $outDir/$sampleID.somaticCNA.varscan.copynumber.called \
#		--output-homdel-file $outDir/$sampleID.somaticCNA.varscan.copynumber.called.homdel \
#		--min-coverage 20 \
#		--min-tumor-coverage 10 \
#		--max-homdel-coverage 5 \
#		--amp-threshold 0.25 \
#		--del-threshold 0.25 \
#		--min-region-size 10
#fi
#
#Rscript $shDir/plot_varScan.R $outDir/$sampleID.somaticCNA.varscan.copynumber.called $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment $sampleID
#awk -F"\t" -v OFS="\t"  '{print $1,$0}' $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment > $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.tmp
#perl $varScanDIR/mergeSegments.pl $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.tmp \
#	--amp-threshold 0.3 \
#	--del-threshold -0.3 \
#	--size-threshold 0.25 \
#	--ref-arm-sizes $varScanDIR/ref-arm-sizes \
#	--output-basename $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge
#rm $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.tmp
#
##filter: 8 continuous bin
#awk -F"\t" -v OFS="\t" 'NR==1 || ($8!="neutral" && $6>=8)' $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.tsv > $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.filter.tsv
#awk -F"\t" -v OFS="\t" 'NR>1{id++;print $1,"varscan",$8,$2,$3,".",".",".","seg_mean="$4";num_segments="$5";num_markers="$6";event_size="$9";size_class="$10";chrom_arm="$11";arm_fraction="$12";chrom_fraction="$13";SVID="id";SVType="$8  }' $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.filter.tsv > $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.filter.gff
#
## annotation
#var_annotation.sh -m "varcan.cnv" $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.filter.gff $sampleID
#
#perl $shDir/somatic_cnv_varScan.stat.pl $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.filter.gff.ann.xls > $outDir/$sampleID.varcan.cnv.ann.summary.xls
#
#awk -F"\t" 'NR==1 || ($27=="deletion"&&($30=="TSG" || $32~/deletion/)) || ($27=="amplification"&&($30=="Oncogene" || $32~/Amplification/))' $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.filter.gff.ann.geneInfo.xls > $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.filter.gff.ann.geneInfo.knownDriver.xls

Rscript $shDir/plot_varScan_23.4.R -b $outDir/$sampleID.somaticCNA.varscan.copynumber.called -s $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment -c $outDir/$sampleID.somaticCNA.varscan.copynumber.called.segment.merge.events.filter.tsv -p $outDir/$sampleID.logR.png

echo end at: `date`
