#!/bin/bash

echo "*** somatic SNV by varScan ***"

TR=""
optTR=""
iniFile="`dirname $0`/../parameter/init_human.sh"
#_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/bwa_0.7.12/human_g1k_v37_decoy.fasta"
tumorFreq=0.1
normalFreq=0.05

while getopts c:r:f:a:b: opt
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
			TR="$OPTARG"
			optTR="-r $TR"
		fi
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
	a)
		normalFreq=$OPTARG
		;;
	b)
		tumorFreq=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-f reference] [-a normal max freq, default 0.05] [-b tumor min freq, default 0.1] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-f reference] [-a normal max freq, default 0.05] [-b tumor min freq, default 0.1] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile
#refData=$_refData
#REF=$_refData

module load samtools/0.1.18

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir

### call ###
if [ ! -f "$outDir/$sampleID.varScan.snp.Somatic.hc.vcf" ]; then
	samtools mpileup -q 1 $optTR -f $refData $normalBam $tumorBam \
		| java -Xmx4G -jar $varScanDIR/VarScan.jar somatic --mpileup 1 \
		--min-coverage-normal 8 \
		--min-coverage-tumor 8 \
		--p-value 0.1 \
		--min-var-freq $tumorFreq \
		--output-snp $outDir/$sampleID.varScan.snp \
		--output-indel $outDir/$sampleID.varScan.indel \
		--output-vcf
	
	java -Xmx4G -jar $varScanDIR/VarScan.jar processSomatic \
		$outDir/$sampleID.varScan.snp.vcf \
		--min-tumor-freq $tumorFreq \
		--max-normal-freq $normalFreq 
	
	java -Xmx4G -jar $varScanDIR/VarScan.jar processSomatic \
		$outDir/$sampleID.varScan.indel.vcf \
		--min-tumor-freq $tumorFreq \
		--max-normal-freq $normalFreq
fi

### fpfilter

## somatic
awk '!/^#/ {print $1,$2,$2}' $outDir/$sampleID.varScan.snp.Somatic.hc.vcf > $outDir/$sampleID.varScan.snp.Somatic.hc.list
bam-readcount $tumorBam -w 3 -q 1 -b 13 -f $refData -l $outDir/$sampleID.varScan.snp.Somatic.hc.list > $outDir/$sampleID.varScan.snp.Somatic.hc.readcount
$varScanDIR/fpfilter_vcf.pl $outDir/$sampleID.varScan.snp.Somatic.hc.vcf $outDir/$sampleID.varScan.snp.Somatic.hc.readcount --output-basename $outDir/$sampleID.varScan.snp.Somatic.hc.fpfilter
awk '/^#/' $outDir/$sampleID.varScan.snp.Somatic.hc.vcf > $outDir/$sampleID.varScan.call.Somatic.vcf
cat $outDir/$sampleID.varScan.snp.Somatic.hc.fpfilter.pass >> $outDir/$sampleID.varScan.call.Somatic.vcf

## germline
awk '!/^#/ {print $1,$2,$2}' $outDir/$sampleID.varScan.snp.Germline.hc.vcf > $outDir/$sampleID.varScan.snp.Germline.hc.list
bam-readcount $tumorBam -w 3 -q 1 -b 13 -f $refData -l $outDir/$sampleID.varScan.snp.Germline.hc.list > $outDir/$sampleID.varScan.snp.Germline.hc.readcount
$varScanDIR/fpfilter_vcf.pl $outDir/$sampleID.varScan.snp.Germline.hc.vcf $outDir/$sampleID.varScan.snp.Germline.hc.readcount --output-basename $outDir/$sampleID.varScan.snp.Germline.hc.fpfilter
awk '/^#/' $outDir/$sampleID.varScan.snp.Germline.hc.vcf > $outDir/$sampleID.varScan.call.Germline.vcf
cat $outDir/$sampleID.varScan.snp.Germline.hc.fpfilter.pass >> $outDir/$sampleID.varScan.call.Germline.vcf

## LOH
awk '!/^#/ {print $1,$2,$2}' $outDir/$sampleID.varScan.snp.LOH.hc.vcf > $outDir/$sampleID.varScan.snp.LOH.hc.list
bam-readcount $tumorBam -w 3 -q 1 -b 13 -f $refData -l $outDir/$sampleID.varScan.snp.LOH.hc.list > $outDir/$sampleID.varScan.snp.LOH.hc.readcount
$varScanDIR/fpfilter_vcf.pl $outDir/$sampleID.varScan.snp.LOH.hc.vcf $outDir/$sampleID.varScan.snp.LOH.hc.readcount --output-basename $outDir/$sampleID.varScan.snp.LOH.hc.fpfilter
awk '/^#/' $outDir/$sampleID.varScan.snp.LOH.hc.vcf > $outDir/$sampleID.varScan.call.LOH.vcf
cat $outDir/$sampleID.varScan.snp.LOH.hc.fpfilter.pass >> $outDir/$sampleID.varScan.call.LOH.vcf

echo end at: `date`
echo "*** Finished somatic SNV by varScan ***"
