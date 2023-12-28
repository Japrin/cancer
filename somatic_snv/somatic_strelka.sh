#!/bin/bash

echo "*** somatic mutation by strelka ***"

TR=""
optTR=""
iniFile="`dirname $0`/../parameter/init_human.hg38.sh"
source $iniFile
_refData=$refData
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

#source $iniFile
refData=$_refData
REF=$_refData

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

##STRELKA_1_0_14_PATH="/WPSnew/zhenglt/01.bin/var/strelka/strelka_workflow-1.0.14/mybuild"

##mkdir -p $outDir
config_file=""

if [ -f "$TR" ]; then
	config_file="$STRELKA_1_0_14_PATH/etc/strelka_config_bwa_default_TR.ini"
else
	config_file="$STRELKA_1_0_14_PATH/etc/strelka_config_bwa_default_WGS.ini"
fi
### call ###
if [ ! -f "$outDir/results/passed.somatic.snvs.vcf" ]; then

	echo "good"
	$STRELKA_1_0_14_PATH/bin/configureStrelkaWorkflow.pl \
		--tumor=$tumorBam \
		--normal=$normalBam \
		--ref=$refData \
		--config=$config_file \
		--output-dir=$outDir
	
	make -C $outDir -j 8

fi

### invcf (bgziped) REF
addMutFreq()
{
    local invcf=$1
    #local REF=$2
    ##bgzip -cd $invcf 
    cat $invcf \
    | perl -ane 'chomp;
                if(/^#/){ print "$_\n"; next; }
                my $line=$_;
                @F=split /\t/; $F[0]=~s/^chr//;
                if($F[0]=~/^([0-9]+|X|Y)$/ && $F[6]=~/PASS/){ print "$line\n"; } ' \
    | bcftools norm -m-both \
    | bcftools norm -f $REF \
    | $shDir/strelka.addFreq.pl
    
    ## ${invcf/.vcf/.flt.vcf}
}


if [ -f "$TR" ]; then
	(
		awk '/^#/' $outDir/results/passed.somatic.snvs.vcf
		bedtools.static.binary intersect -a $outDir/results/passed.somatic.snvs.vcf   -b $TR | uniq 
	) | addMutFreq > $outDir/results/$sampleID.strelka.somatic.snvs.filter.vcf

	(
		awk '/^#/' $outDir/results/passed.somatic.indels.vcf
		bedtools.static.binary intersect -a $outDir/results/passed.somatic.indels.vcf -b $TR | uniq 
	) | addMutFreq > $outDir/results/$sampleID.strelka.somatic.indels.filter.vcf
else
	cat $outDir/results/passed.somatic.snvs.vcf   | addMutFreq > $outDir/results/$sampleID.strelka.somatic.snvs.filter.vcf
	cat $outDir/results/passed.somatic.indels.vcf | addMutFreq > $outDir/results/$sampleID.strelka.somatic.indels.filter.vcf
fi

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

run_annovar $outDir/results/$sampleID.strelka.somatic.snvs.filter.vcf "hg38" $HumanDB
run_annovar $outDir/results/$sampleID.strelka.somatic.indels.filter.vcf "hg38" $HumanDB



#####if [ "$optZ" = "Y" ] || [ "$optZ" = "y" ]
#####then
#####	var_annotation.sh -c $iniFile -a somatic -b SNP   -m strelka.snv   $optG $outDir/results/$sampleID.strelka.somatic.snvs.filter.vcf    N,T
#####	var_annotation.sh -c $iniFile -a somatic -b INDEL -m strelka.indel $optG $outDir/results/$sampleID.strelka.somatic.indels.filter.vcf  N,T
#####
#####	### to maf format
#####	bgzip -cd $outDir/results/$sampleID.strelka.somatic.snvs.filter.reformated.vcf.gz \
#####		| python2 $PIPELINE/cancer/var/vcf2mafv1.31.strelka.py -t - -s "WES" -a $sampleID \
#####			-g $HumanDB/hg19_knownGene.SymbolToLocus.txt \
#####			-o $outDir/results/$sampleID.strelka.somatic.snvs.filter.reformated.maf
#####
#####	bgzip -cd $outDir/results/$sampleID.strelka.somatic.indels.filter.reformated.vcf.gz \
#####		| python2 $PIPELINE/cancer/var/vcf2mafv1.31.strelka.py -t - -s "WES" -a $sampleID \
#####			-g $HumanDB/hg19_knownGene.SymbolToLocus.txt \
#####			-o $outDir/results/$sampleID.strelka.somatic.indels.filter.reformated.maf
#####fi

echo end at: `date`
echo "*** Finished somatic mutation by strelka ***"
