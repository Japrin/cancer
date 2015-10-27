#!/bin/bash -eu
echo "*** variant calling by GATK ***"

TR=""
optTR=""
iniFile="/WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/parameter/init_human.sh"

opt_P=""

while getopts c:r:p: opt
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
			optTR="-L $TR"
		else
			echo "WARNING: invalid target file ($OPTARG), no target will be used"
		fi
		;;
	p)
		opt_P="$OPTARG";
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-p extra UG parameters ] <sampleID> <inbam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-p extra UG parameters ] <sampleID> <inbam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile 

sampleID=$1
in_bam=$2
outDir=$3

mkdir -p $outDir

raw_vcf=$outDir/$sampleID.GATK.var.raw.vcf
flt_vcf=$outDir/$sampleID.GATK.var.flt.vcf
snp_vcf=$outDir/$sampleID.GATK.snp.vcf
snp_flt_vcf=$outDir/$sampleID.GATK.snp.flt.vcf
indel_vcf=$outDir/$sampleID.GATK.indel.vcf
indel_flt_vcf=$outDir/$sampleID.GATK.indel.flt.vcf

if [ -f $snp_flt_vcf ]
then
    echo "## $snp_flt_vcf exists, skip this step"
	exit 0
fi

echo ">>> Running the unified genotyper for SNP/INDEL calling"
echo ">>> INDEL"
java -Xms5g -Xmx5g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	-T UnifiedGenotyper \
	-I $in_bam \
	-R $REF \
	-o $indel_vcf $optTR \
	-dcov 1000 \
	-A QualByDepth \
	-A FisherStrand \
	-A AlleleBalance \
	-A Coverage \
	-A MappingQualityZero \
	-A TandemRepeatAnnotator \
	-A VariantType \
	-A DepthPerAlleleBySample \
	-baq CALCULATE_AS_NECESSARY \
 	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0 \
	-glm INDEL \
	-et NO_ET \
	-rf BadCigar \
	-K $GATKKey $opt_P

filter="QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0"

java -Xms5g -Xmx5g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $REF \
	-o $indel_flt_vcf \
	-V $indel_vcf \
	--filterExpression "$filter" \
	--filterName "StandardFilter" \
	-et NO_ET \
	-rf BadCigar \
	-K $GATKKey
#perl -i -F"\t" -ane 'if(/^#/||$F[6] eq "."||$F[6] eq "PASS"){print}' $indel_flt_vcf
rm $indel_flt_vcf.idx

echo ">>> SNP"
java -Xms5g -Xmx5g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	-T UnifiedGenotyper \
	-I $in_bam \
	-R $REF \
	-o $snp_vcf $optTR \
	-dcov 1000 \
	-A QualByDepth \
	-A FisherStrand \
	-A AlleleBalance \
	-A Coverage \
	-A MappingQualityZero \
	-A TandemRepeatAnnotator \
	-A VariantType \
	-A DepthPerAlleleBySample \
	-baq CALCULATE_AS_NECESSARY \
 	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0 \
	-glm SNP \
	-et NO_ET \
	-rf BadCigar \
	-K $GATKKey $opt_P

filter="QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
###filter="QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

java -Xms5g -Xmx5g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $REF \
	-o $snp_flt_vcf \
	-V $snp_vcf \
	--filterExpression "$filter" \
	--filterName "StandardFilter" \
	-et NO_ET \
	-rf BadCigar \
	-K $GATKKey
#perl -i -F"\t" -ane 'if(/^#/||$F[6] eq "."||$F[6] eq "PASS"){print}' $snp_flt_vcf
rm $snp_flt_vcf.idx

echo ""

echo "*** Finished SNP Analysis using the GATK Unified Genotyper ***"
