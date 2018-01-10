#!/bin/bash

echo "*** somatic mutation by strelka2 ***"

TR=""
optTR=""
iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
#_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optZ="N"
optG=""

shDir=`dirname $0`

while getopts c:r:f:z:g: opt
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
	g)
		if [ -f $OPTARG ];then
			optG=" -g $OPTARG "
		fi
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
	z)
		optZ="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-g gender file] [-f reference] [-z N|Y] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-g gender file] [-f reference] [-z N|Y] <sampleID> <normalBam> <tumorBam> <outDir>"
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

mkdir -p $outDir/manta
mkdir -p $outDir/strelka
mkdir -p $outDir/annovar

#mkdir -p $outDir
config_file=""

STRELKA_INSTALL_PATH="/Share/BP/zhenglt/01.bin/strelka/strelka-2.8.3.centos5_x86_64"
MANTA_INSTALL_PATH="/Share/BP/zhenglt/01.bin/strelka/manta-1.2.1.centos6_x86_64"

####### Manta run #######
### --outputContig
if [ -f "$TR" ]; then
    ${MANTA_INSTALL_PATH}/bin/configManta.py \
        --normalBam $normalBam \
        --tumorBam $tumorBam \
        --referenceFasta $REF \
        --runDir $outDir/manta \
        --exome \
        --callRegions $TR
else
    ${MANTA_INSTALL_PATH}/bin/configManta.py \
        --normalBam $normalBam \
        --tumorBam $tumorBam \
        --referenceFasta $REF \
        --runDir $outDir/manta
fi
$outDir/manta/runWorkflow.py -m local -j 8

if [ -f "$TR" ]; then
    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
        --normalBam $normalBam \
        --tumorBam $tumorBam \
        --referenceFasta $REF \
        --indelCandidates $outDir/manta/results/variants/candidateSmallIndels.vcf.gz \
        --runDir $outDir/strelka \
        --outputCallableRegions \
        --exome \
        --callRegions $TR
else
    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
        --normalBam $normalBam \
        --tumorBam $tumorBam \
        --referenceFasta $REF \
        --indelCandidates $outDir/manta/results/variants/candidateSmallIndels.vcf.gz \
        --runDir $outDir/strelka \
        --outputCallableRegions
fi
$outDir/strelka/runWorkflow.py -m local -j 8

###### add mutation frequency

bgzip -cd $outDir/strelka/results/variants/somatic.snvs.vcf.gz \
    | perl -ane 'chomp;
                if(/^#/){ print "$_\n"; next; }
                @F=split /\t/; $F[0]=~s/^chr//;
                if($F[0]=~/^([0-9]+|X|Y)$/ && $F[6]=~/PASS/){ print "$_\n"; } ' \
    | $shDir/strelka.addFreq.pl \
    > $outDir/annovar/$sampleID.strelka.somatic.snvs.filter.vcf
bgzip -cd $outDir/strelka/results/variants/somatic.indels.vcf.gz \
    | perl -ane 'chomp;
                if(/^#/){ print "$_\n"; next; }
                @F=split /\t/; $F[0]=~s/^chr//;
                if($F[0]=~/^([0-9]+|X|Y)$/ && $F[6]=~/PASS/){ print "$_\n"; } ' \
    | $shDir/strelka.addFreq.pl \
    > $outDir/annovar/$sampleID.strelka.somatic.indels.filter.vcf

if [ "$optZ" = "Y" ] || [ "$optZ" = "y" ]
then
	var_annotation.sh -c $iniFile -a somatic -b SNP   -m strelka.snv   $optG $outDir/annovar/$sampleID.strelka.somatic.snvs.filter.vcf    N,T
	var_annotation.sh -c $iniFile -a somatic -b INDEL -m strelka.indel $optG $outDir/annovar/$sampleID.strelka.somatic.indels.filter.vcf  N,T
	### to maf format
	bgzip -cd $outDir/annovar/$sampleID.strelka.somatic.snvs.filter.reformated.vcf.gz \
		| $PIPELINE/cancer/var/vcf2mafv1.31.strelka.py -t - -s "WES" -a $sampleID \
			-g $HumanDB/hg19_knownGene.SymbolToLocus.txt \
			-o $outDir/annovar/$sampleID.strelka.somatic.snvs.filter.reformated.maf
	bgzip -cd $outDir/annovar/$sampleID.strelka.somatic.indels.filter.reformated.vcf.gz \
		| $PIPELINE/cancer/var/vcf2mafv1.31.strelka.py -t - -s "WES" -a $sampleID \
			-g $HumanDB/hg19_knownGene.SymbolToLocus.txt \
			-o $outDir/annovar/$sampleID.strelka.somatic.indels.filter.reformated.maf
fi

echo end at: `date`
echo "*** Finished somatic mutation by strelka ***"
