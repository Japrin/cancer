#!/bin/bash

echo "*** somatic mutation by strelka2 ***"

TR=""
optTR=""
iniFile="`dirname $0`/../parameter/init_human.hg38.sh"
source $iniFile
_refData=$refData
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

##source $iniFile
refData=$_refData
REF=$_refData

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir/manta
mkdir -p $outDir/strelka
mkdir -p $outDir/annovar

#mkdir -p $outDir
config_file=""

####### Manta run #######

if [ ! -f "$outDir/manta/results/variants/candidateSmallIndels.vcf.gz" ]; then
    if [ -f "$TR" ]; then
        ${MANTA_INSTALL_PATH}/bin/configManta.py \
            --normalBam $normalBam \
            --tumorBam $tumorBam \
            --referenceFasta $REF \
            --runDir $outDir/manta \
            --exome \
            --generateEvidenceBam \
            --callRegions $TR
    else
        ${MANTA_INSTALL_PATH}/bin/configManta.py \
            --normalBam $normalBam \
            --tumorBam $tumorBam \
            --referenceFasta $REF \
            --runDir $outDir/manta
    fi
    $outDir/manta/runWorkflow.py -m local -j 8
fi

if [ ! -f "$outDir/strelka/results/variants/somatic.snvs.vcf.gz" ]; then
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
fi

###### add mutation frequency

### invcf (bgziped) REF
addMutFreq()
{
    local invcf=$1
    local REF=$2
    bgzip -cd $invcf \
    | perl -ane 'chomp;
                if(/^#/){ print "$_\n"; next; }
                my $line=$_;
                @F=split /\t/; $F[0]=~s/^chr//;
                if($F[0]=~/^([0-9]+|X|Y)$/ && $F[6]=~/PASS/){ print "$line\n"; } ' \
    | bcftools norm -m-both \
    | bcftools norm -f $REF \
    | $shDir/strelka.addFreq.pl \
    > ${invcf/.vcf.gz/.flt.vcf}
}

addMutFreq $outDir/strelka/results/variants/somatic.snvs.vcf.gz $refData
addMutFreq $outDir/strelka/results/variants/somatic.indels.vcf.gz $refData

### annovarFile  annVer HumanDB
run_annovar()
{
    local annovarFile=$1
    local annVer=$2
    local HumanDB=$3

    local nvariant=$(awk '!/^#/' $annovarFile | wc -l)

    if [ "$nvariant" -gt 0 ];then

        table_annovar.pl $annovarFile $HumanDB -buildver $annVer -otherinfo -remove -nastring . \
                                -protocol ensGene,cytoBand,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp150,ljb26_all,cosmic97 \
                                -operation g,r,f,f,f,f,f,f,f,f \
                                -vcfinput
    
        ## filter by population
        perl -F"\t" -ane 'if(/^#/ || ( ($F[7]!~/snp150=rs/ && $F[2] eq".") ||  ($F[7]=~/cosmic\d\d=ID/ && ($F[7]=~/1000g2015aug_all=\./ || ($F[7]=~/1000g2015aug_all=(.+?);/ && $1<0.01) ) ) ) ){print}' \
            $annovarFile.hg38_multianno.vcf \
            > $annovarFile.hg38_multianno.flt.vcf
        bgzip -f $annovarFile.hg38_multianno.flt.vcf
        tabix -f -p vcf $annovarFile.hg38_multianno.flt.vcf.gz
    
        perl -F"\t" -ane 'if(/^Chr/ || ( ($F[16]!~/^rs/) || ($F[42]=~/^ID/ && ($F[12] eq "." || $F[12]<0.01) ) ) ){print}' \
            $annovarFile.hg38_multianno.txt \
            > $annovarFile.hg38_multianno.flt.txt
    else
        echo "Only $nvariant variansts detected! Annotation will not run! "
    fi

}

run_annovar $outDir/strelka/results/variants/somatic.snvs.flt.vcf "hg38" $HumanDB
run_annovar $outDir/strelka/results/variants/somatic.indels.flt.vcf "hg38" $HumanDB


#####if [ "$optZ" = "Y" ] || [ "$optZ" = "y" ]
#####then
#####	var_annotation.sh -c $iniFile -a somatic -b SNP   -m strelka.snv   $optG $outDir/annovar/$sampleID.strelka.somatic.snvs.filter.vcf    N,T
#####	var_annotation.sh -c $iniFile -a somatic -b INDEL -m strelka.indel $optG $outDir/annovar/$sampleID.strelka.somatic.indels.filter.vcf  N,T
#####	### to maf format
#####	bgzip -cd $outDir/annovar/$sampleID.strelka.somatic.snvs.filter.reformated.vcf.gz \
#####		| $PIPELINE/cancer/var/vcf2mafv1.31.strelka.py -t - -s "WES" -a $sampleID \
#####			-g $HumanDB/hg19_knownGene.SymbolToLocus.txt \
#####			-o $outDir/annovar/$sampleID.strelka.somatic.snvs.filter.reformated.maf
#####	bgzip -cd $outDir/annovar/$sampleID.strelka.somatic.indels.filter.reformated.vcf.gz \
#####		| $PIPELINE/cancer/var/vcf2mafv1.31.strelka.py -t - -s "WES" -a $sampleID \
#####			-g $HumanDB/hg19_knownGene.SymbolToLocus.txt \
#####			-o $outDir/annovar/$sampleID.strelka.somatic.indels.filter.reformated.maf
#####fi

echo end at: `date`
echo "*** Finished somatic mutation by strelka ***"
