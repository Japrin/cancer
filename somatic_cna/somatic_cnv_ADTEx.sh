#!/bin/bash

iniFile="`dirname $0`/../parameter/init_human.sh"
#_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
TR=""
normalSNPFile=""
gender="M"
opt_e="--estimatePloidy"

while getopts c:b:r:f:g:e: opt
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
			echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] [-f reference] [-g gender,default \"M\"] [-e ADTex option, default \"--estimatePloidy\" ] [-b normal snp file, called by GATK] <sampleID> <normalBam> <tumorBam> <outDir>"
			exit 1
		fi
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
	b)
		normalSNPFile=$OPTARG
		;;
	g)
		gender=$OPTARG
		;;
	e)
		opt_e=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] [-f reference] [-g gender,default \"M\"] [-e ADTex option, default \"--estimatePloidy\" ] [-b normal snp file, called by GATK] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion (required)] [-f reference] [-g gender,default \"M\"] [-e ADTex option, default \"--estimatePloidy\" ] [-b normal snp file, called by GATK] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile
#refData=$_refData
#REF=$_refData

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

shDir=`dirname $0`
mkdir -p $outDir

#ADTExDir=/Share/BP/zhenglt/01.bin/ADTEx/ADTEx.v.1.0.4
ADTExDir=/WPSnew/zhenglt/01.bin/var/adtex/ADTEx.v.1.0.4
module load samtools/0.1.19
## bedtools v2.24 and above are not compatible with ADTEx now
#module unload bedtools/v2.25.0
module unload bedtools/2.27.1
module load bedtools/2.23.0
export TMPDIR="$outDir"

### get the BAF file
if [ ! -f $outDir/$sampleID.hetSNV.baf ];then
	if [ -f $normalSNPFile ];then
		bgzip -cd $normalSNPFile \
			| $PIPELINE/cancer/rna/mbased.filterSNV.pl -n -f -a 5 -b 0.1 -g $gender \
			> $outDir/$sampleID.hetSNV.bed
		samtools mpileup -f $refData -l $outDir/$sampleID.hetSNV.bed $normalBam $tumorBam > $outDir/$sampleID.pileup
		countMultiPileup.pl -maxNormal 1 -minDepth 1 $outDir/$sampleID.pileup > $outDir/$sampleID.hetSNV.RC
		awk -F"\t" -v OFS="\t" 'BEGIN{print "chrom\tSNP_loc\tcontrol_BAF\ttumor_BAF\tmirrored_BAF\tcontrol_doc\ttumor_doc\tcontrol_vdoc\ttumor_vdoc"} $6+$7>=10{ t_BAF=$7/($7+$6);print $1,$2,$5/($4+$5),t_BAF,(t_BAF>0.5?(1-t_BAF):t_BAF),$4+$5,$6+$7,$5,$7 }' $outDir/$sampleID.hetSNV.RC > $outDir/$sampleID.hetSNV.baf
	fi
fi

if [ ! -f $outDir/$sampleID.tumor.coverage.gz ];then
	## get the DOC file
	samtools view -uF 0x400 $normalBam | coverageBed -abam - -b $TR -d | gzip -c > $outDir/$sampleID.normal.coverage.gz
	samtools view -uF 0x400 $tumorBam  | coverageBed -abam - -b $TR -d | gzip -c > $outDir/$sampleID.tumor.coverage.gz
fi

## run ADTEx
if [ ! -f $outDir/out/cnv.result ];then
	####rm -r $outDir/out
	if [ -f $outDir/$sampleID.hetSNV.baf ];then
		python $ADTExDir/ADTEx.py \
			-n $outDir/$sampleID.normal.coverage.gz \
			-t $outDir/$sampleID.tumor.coverage.gz \
			-b $TR \
			-o $outDir/out \
			--baf $outDir/$sampleID.hetSNV.baf \
			$opt_e --DOC
			#$opt_e -p --DOC
	else
	
		python $ADTExDir/ADTEx.py \
			-n $outDir/$sampleID.normal.coverage.gz \
			-t $outDir/$sampleID.tumor.coverage.gz \
			-b $TR \
			-o $outDir/out \
			$opt_e --DOC
			#$opt_e -p --DOC
	fi
fi

awk -F"\t" -v OFS="\t" 'NR>1{print $1,$8,$9,$10,$4 }' $outDir/out/cnv.result | uniq -c \
	| awk  -v OFS="\t" 'BEGIN{print "chr\tbeg\tend\tratio\tcopy_number\ttarget_number"} {print $2,$3,$4,$5,$6,$1}' \
	> $outDir/out/cnv.result.segment
## ploidy
awk  'BEGIN {l=0;s=0} NR>1{ l+=$3-$2;s+=($3-$2)*$5; } END{print "'$sampleID'\t"l"\t"(s/l)}' $outDir/out/cnv.result.segment > $outDir/out/cnv.ploidy
pl=`cut -f 3 $outDir/out/cnv.ploidy`
## select for those CN != ploidy
awk 'NR==1 || ($5!=int('$pl'+0.5) && $6>=10)' $outDir/out/cnv.result.segment > $outDir/out/cnv.result.segment.CNV.txt

## annotation
awk -F"\t" -v OFS="\t" 'NR>1{ svtype=($4<1?"del":"amp"); print $1,"ADTEx",svtype,$2,$3,".",".",".","ratio="$4";CN="$5";target_number="$6";SVID='$sampleID'.CNV"NR";SVType="svtype";ploidy='$pl'"; }' \
	$outDir/out/cnv.result.segment.CNV.txt \
	> $outDir/out/cnv.result.segment.CNV.gff
var_annotation.sh $outDir/out/cnv.result.segment.CNV.gff $sampleID


if [ -f $outDir/$sampleID.hetSNV.baf ];then
	## LOH
	$shDir/ADTEx.zygosity.segments.pl $outDir/out/zygosity/zygosity.res > $outDir/out/zygosity/zygosity.res.segment
	awk 'NR==1 || $4=="LOH"' $outDir/out/zygosity/zygosity.res.segment > $outDir/out/zygosity/zygosity.res.segment.LOH.txt
	
	## plot
	printf "execute: \n$PIPELINE/cancer/somatic_cna/plot_ADTEx.R $sampleID $outDir/out/cnv.result $outDir/out/cnv.result.segment $outDir/out $outDir/out/zygosity/zygosity.res $outDir/out/zygosity/zygosity.res.segment\n"
	$PIPELINE/cancer/somatic_cna/plot_ADTEx.R $sampleID $outDir/out/cnv.result $outDir/out/cnv.result.segment $outDir/out $outDir/out/zygosity/zygosity.res $outDir/out/zygosity/zygosity.res.segment
    $PIPELINE/cancer/somatic_cna/plot_ADTEx.tmpPDF.R $sampleID $outDir/out/cnv.result $outDir/out/cnv.result.segment $outDir/out $outDir/out/zygosity/zygosity.res $outDir/out/zygosity/zygosity.res.segment
	
	awk -F"\t" -v OFS="\t" 'NR>1{print $1,"ADTEx","LOH",$2,$3,".\t.\t.","mirror_baf="$5";CN="$6";Z="$7";target_number="$8";SVID='$sampleID'.LOH"NR";SVType=LOH;ploidy='$pl'"; }' \
		$outDir/out/zygosity/zygosity.res.segment.LOH.txt \
		> $outDir/out/zygosity/zygosity.res.segment.LOH.gff
	var_annotation.sh $outDir/out/zygosity/zygosity.res.segment.LOH.gff $sampleID

fi

echo end at: `date`
