#!/bin/bash

echo "*** variant calling by varScan ***"

TR=""
optTR=""
iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/bwa_0.7.12/human_g1k_v37_decoy.fasta"
optG="M"
optA="n"

while getopts c:f:r:a:g: opt
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
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
	a)
		optA=$OPTARG
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
	g)
		optG=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-f reference] [-r targetRegion] [-g gender,F for female] [-a yes/no for annotation] <sampleID> <inbam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-f reference] [-r targetRegion] [-g gender,F for female] [-a y/n for annotation] <sampleID> <inbam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile
refData=$_refData
REF=$_refData

module load samtools/0.1.18


sampleID=$1
in_bam=$2
outDir=$3

mkdir -p $outDir

if [ ! -f "$outDir/$sampleID.var.varScan.snp.fpfilter.pass" ];then
	## snp
	samtools mpileup -q 1 $optTR -B -f $refData $in_bam \
		| java -Xmx4G -jar $varScanDIR/VarScan.jar mpileup2snp --output-vcf  \
			--min-coverage 10 \
			--min-reads2 2 \
		        --min-avg-qual 15 \
			--min-var-freq 0.05 \
			--min-freq-for-hom 0.75 \
			--p-value  0.05 \
			--strand-filter 1 \
			> $outDir/$sampleID.var.varScan.snp.vcf
	
	awk '!/^#/ {print $1,$2,$2}' $outDir/$sampleID.var.varScan.snp.vcf > $outDir/$sampleID.var.varScan.snp.list
	bam-readcount $in_bam -w 50 -q 1 -b 13 -f $refData -l $outDir/$sampleID.var.varScan.snp.list > $outDir/$sampleID.var.varScan.snp.readcount
	
	java -Xmx4G -jar $varScanDIR/VarScan.jar fpfilter \
		$outDir/$sampleID.var.varScan.snp.vcf \
		$outDir/$sampleID.var.varScan.snp.readcount \
		--min-ref-basequal 20 \
		--min-var-basequal 20 \
	       	--output-file $outDir/$sampleID.var.varScan.snp.fpfilter.pass  \
		--filtered-file $outDir/$sampleID.var.varScan.snp.fpfilter.fail
	
fi


if [ ! -f "$outDir/$sampleID.var.varScan.indel.fpfilter.pass" ];then
	## indel
	samtools mpileup -q 1 $optTR -B -f $refData $in_bam \
		| java -Xmx4G -jar $varScanDIR/VarScan.jar mpileup2indel --output-vcf  \
			--min-coverage 10 \
			--min-reads2 2 \
		        --min-avg-qual 15 \
			--min-var-freq 0.05 \
			--min-freq-for-hom 0.75 \
			--p-value  0.05 \
			--strand-filter 1 \
			> $outDir/$sampleID.var.varScan.indel.vcf
	
	awk '!/^#/ {print $1,$2,$2}' $outDir/$sampleID.var.varScan.indel.vcf > $outDir/$sampleID.var.varScan.indel.list
	bam-readcount $in_bam -w 50 -q 1 -b 13 -f $refData -l $outDir/$sampleID.var.varScan.indel.list > $outDir/$sampleID.var.varScan.indel.readcount
	
	ln -s $outDir/$sampleID.var.varScan.indel.vcf $outDir/$sampleID.var.varScan.indel.fpfilter.pass
	#java -Xmx4G -jar $varScanDIR/VarScan.jar fpfilter \
	#	$outDir/$sampleID.var.varScan.indel.vcf \
	#	$outDir/$sampleID.var.varScan.indel.readcount \
	#	--min-ref-basequal 20 \
	#	--min-var-basequal 20 \
	#      	--output-file $outDir/$sampleID.var.varScan.indel.fpfilter.pass  \
	#	--filtered-file $outDir/$sampleID.var.varScan.indel.fpfilter.fail
	
fi

awk ' {  if(/^##/){print $0}else { if(/^#CHROM/) { line_chr=$0; flag=0 }}  if(!/^#/){ if(flag==0){print line_chr;flag=1} if($1!~/^(GL|hs)/){ print $0 }  }  }' $outDir/$sampleID.var.varScan.snp.fpfilter.pass > $outDir/$sampleID.var.varScan.snp.call.vcf
awk ' {  if(/^##/){print $0}else { if(/^#CHROM/) { line_chr=$0; flag=0 }}  if(!/^#/){ if(flag==0){print line_chr;flag=1} if($1!~/^(GL|hs)/){ print $0 }  }  }' $outDir/$sampleID.var.varScan.indel.fpfilter.pass > $outDir/$sampleID.var.varScan.indel.call.vcf

if [ "$optA" == "y" ];then
	var_annotation.sh -m varScan -g $optG $outDir/$sampleID.var.varScan.snp.call.vcf $sampleID
	var_annotation.sh -m varScan -g $optG $outDir/$sampleID.var.varScan.indel.call.vcf $sampleID
fi

echo end at: `date`
