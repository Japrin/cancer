#!/bin/bash

echo "*** variant calling by samtools ***"

TR=""
optTR=""
iniFile="`dirname $0`/../parameter/init_human.sh"
optG="-v"
optA=""

while getopts c:r:a:g opt
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
	a)
		optA="-r $OPTARG"
		;;
	g)
		optG=""
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-a chr] [-g genotyping] <sampleID> <inbam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-a chr] [-g genotyping] <sampleID> <inbam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile 

sampleID=$1
in_bam=$2
outDir=$3

raw_bcf=$outDir/$sampleID.samtools.var.raw.bcf
raw_vcf=$outDir/$sampleID.samtools.var.raw.vcf
flt_vcf=$outDir/$sampleID.samtools.var.flt.vcf
snp_vcf=$outDir/$sampleID.samtools.snp.vcf
indel_vcf=$outDir/$sampleID.samtools.indel.vcf

mkdir -p $outDir

if [ ! -f "$flt_vcf.gz" ]
then

	samtools mpileup $optTR $optA -q 1 -C 50 -m 2 -F 0.002 -t DP -t SP -t DV -t DP4 -go $raw_bcf -f $refData $in_bam && \
	bcftools call -vmO z -o $raw_vcf.gz $raw_bcf
	bcftools filter -s LowQual -e '%QUAL<20 || MQ<10 || DP<4 || DP>5000 || DV<2' --SnpGap 3 --IndelGap 2 -O z -o $flt_vcf.gz $raw_vcf.gz
	tabix -f -p vcf $flt_vcf.gz
	bcftools filter -i 'Type="snp"' -O z -o $snp_vcf.gz $flt_vcf.gz
	bcftools filter -i 'Type="indel"' -O z -o $indel_vcf.gz $flt_vcf.gz
	tabix -f -p vcf $snp_vcf.gz
	tabix -f -p vcf $indel_vcf.gz
	
	#bcftools view $raw_vcf.gz | vcfutils.pl varFilter -d 4 -D 1000000000 - > $flt_vcf
		
	#samtools mpileup $optTR -q 1 -C 50 -D -S -m 2 -F 0.002 -uf $refData $in_bam | bcftools view $optTR -bcg $optG - > $raw_bcf && \
	#bcftools view $raw_bcf | vcfutils.pl varFilter -d 4 -D 1000000000 - > $flt_vcf
	#bgzip -f $flt_vcf && \
	#tabix -f -p vcf $flt_vcf.gz

	#bgzip -cd $flt_vcf.gz | awk '/^#/ || /INDEL;/' | bgzip -c > $indel_vcf.gz && \
	#tabix -f -p vcf $indel_vcf.gz
	#bgzip -cd $flt_vcf.gz | awk '/^#/ || !/INDEL;/' | bgzip -c > $snp_vcf.gz && \
	#tabix -f -p vcf $snp_vcf.gz

fi

if [ "$optG" = "-v" ]
then
	echo "variant only"
else
	echo "genotype calling"
	bcftools view $raw_bcf | awk '/^#/||$6>=20' | sed 's/;\+/;/' > $outDir/$sampleID.samtools.hc.genotype.vcf && \
	bgzip -f $outDir/$sampleID.samtools.hc.genotype.vcf && \
	tabix -f -p vcf $outDir/$sampleID.samtools.hc.genotype.vcf.gz
fi

echo finish variant calling at: `date`

echo end at: `date`
echo "*** Finished variant calling by samtools ***"
