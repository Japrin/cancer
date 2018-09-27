#!/bin/bash -eu
echo "*** variant calling by GATK ***"

TR=""
optTR=""
iniFile="`dirname $0`/../parameter/init_human.sh"

optM=6
opt_P=""

while getopts c:r:p:m: opt
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
		TR="$OPTARG"
		optTR="-L $TR"
		;;
	m)
		optM=$OPTARG
		;;
	p)
		opt_P="$OPTARG";
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-p extra HC parameters ] [-m memrory(GB), default 6] <sampleID> <inbam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-p extra HC parameters ] [-m memrory(GB), default 6] <sampleID> <inbam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile 

sampleID=$1
in_bam=$2
outDir=$3

raw_vcf=$outDir/$sampleID.GATK.var.raw.vcf
flt_vcf=$outDir/$sampleID.GATK.var.flt.vcf
t_vcf=$outDir/$sampleID.GATK.var.t.vcf
recal_SNP_file=$outDir/$sampleID.GATK.SNP.recal
tranches_SNP_File=$outDir/$sampleID.GATK.SNP.tranches
recal_INDEL_file=$outDir/$sampleID.GATK.INDEL.recal
tranches_INDEL_File=$outDir/$sampleID.GATK.INDEL.tranches

snp_vcf=$outDir/$sampleID.GATK.snp.vcf
snp_flt_vcf=$outDir/$sampleID.GATK.snp.flt.vcf
indel_vcf=$outDir/$sampleID.GATK.indel.vcf
indel_flt_vcf=$outDir/$sampleID.GATK.indel.flt.vcf


echo ">>> Running the haplotype caller for SNP/INDEL calling"

if [ ! -f "$raw_vcf" ];then
	java -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R $REF \
		-I $in_bam \
		--dbsnp $knownSites \
		-stand_call_conf 30.0 \
		-stand_emit_conf 10.0 \
		--genotyping_mode DISCOVERY \
		-o $raw_vcf $optTR $opt_P 
		##-et NO_ET \
		##-K $GATKKey
fi


java -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	-T VariantRecalibrator \
	-R $REF \
	-input $raw_vcf \
	-recalFile $recal_SNP_file \
	--maxGaussians 4 --minNumBadVariants 3000 \
	-tranchesFile $tranches_SNP_File \
	-nt 4 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $bundleDir/hapmap_3.3.b37.vcf \
	-resource:omni,known=false,training=true,truth=true,prior=12.0  $bundleDir/1000G_omni2.5.b37.vcf \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $bundleDir/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $bundleDir/dbsnp_137.b37.vcf \
	-an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP \
	-mode SNP
	##-et NO_ET \
	##-K $GATKKey

java -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	-T VariantRecalibrator \
	-R $REF \
	-input $raw_vcf \
	-recalFile $recal_INDEL_file \
	--maxGaussians 4 --minNumBadVariants 3000 \
	-tranchesFile $tranches_INDEL_File \
	-nt 4 \
	-resource:mills,known=false,training=true,truth=true,prior=12.0 $bundleDir/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $bundleDir/dbsnp_137.b37.vcf \
	-an DP -an FS -an ReadPosRankSum -an MQRankSum \
	-mode INDEL
	##-et NO_ET \
	##-K $GATKKey

java -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	-T ApplyRecalibration \
	-R $REF \
	-input $raw_vcf \
	-recalFile $recal_SNP_file \
	-tranchesFile $tranches_SNP_File \
	-o $t_vcf \
	--ts_filter_level 99.9 \
	-mode SNP
	##-et NO_ET \
	##-K $GATKKey

java -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	-T ApplyRecalibration \
	-R $REF \
	-input $t_vcf \
	-recalFile $recal_INDEL_file \
	-tranchesFile $tranches_INDEL_File \
	-o $flt_vcf \
	--ts_filter_level 99.9 \
	-mode INDEL
	##-et NO_ET \
	##-K $GATKKey

echo end at: `date`

echo "*** Finished SNP Analysis using the GATK Unified Genotyper ***"
