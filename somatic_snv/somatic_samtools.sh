#!/bin/bash

echo "*** somatic SNV by samtools ***"

TR=""
optTR=""
iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
optZ="N"

while getopts c:r:z: opt
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
	z)
		optZ="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [z annotation, Y or N] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [z annotation, Y or N] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir

module load samtools/0.1.19

## raw call
if [ ! -f "$outDir/$sampleID.samtools.somatic.snv.fpfilter.vcf" ]; then
	samtools mpileup $optTR -q 1 -C 50 -m 2 -F 0.002  -gDSf $refData $normalBam $tumorBam > $outDir/$sampleID.mpileup && \
	bcftools view -vcgT pair $outDir/$sampleID.mpileup > $outDir/$sampleID.mpileup.vcf
	#bcftools call -vmO v -o $outDir/$sampleID.mpileup.vcf $outDir/$sampleID.mpileup
	
	## filter
	filter_somatic_snv.pl -normalDP 8 -tumorDP 8 \
	  $outDir/$sampleID.mpileup.vcf > $outDir/$sampleID.samtools.somatic.vcf
	
	awk '/^#/ || /INDEL;/' $outDir/$sampleID.samtools.somatic.vcf   > $outDir/$sampleID.samtools.somatic.indel.vcf
	awk '/^#/ || !/INDEL;/' $outDir/$sampleID.samtools.somatic.vcf  > $outDir/$sampleID.samtools.somatic.snv.vcf
	
	## fpfilter, same as varScan ##
	
	awk '!/^#/ {print $1,$2,$2}' $outDir/$sampleID.samtools.somatic.snv.vcf > $outDir/$sampleID.samtools.somatic.snv.vcf.list
	bam-readcount $tumorBam -q 1 -b 13 -f $refData -l $outDir/$sampleID.samtools.somatic.snv.vcf.list > $outDir/$sampleID.samtools.somatic.snv.vcf.readcount
	$varScanDIR/fpfilter_vcf.pl $outDir/$sampleID.samtools.somatic.snv.vcf $outDir/$sampleID.samtools.somatic.snv.vcf.readcount --output-basename $outDir/$sampleID.samtools.somatic.snv.vcf.fpfilter
	ln -f $outDir/$sampleID.samtools.somatic.snv.vcf.fpfilter.pass $outDir/$sampleID.samtools.somatic.snv.fpfilter.vcf
fi

if [ "$optZ" = "Y" ] || [ "$optZ" = "y" ]
then
	var_annotation.sh -c $iniFile -a somatic -b SNP $outDir/$sampleID.samtools.somatic.snv.fpfilter.vcf N,T
	var_annotation.sh -c $iniFile -a somatic -b INDEL $outDir/$sampleID.samtools.somatic.indel.vcf N,T
fi

echo end at: `date`
echo "*** Finished somatic SNV by samtools ***"
