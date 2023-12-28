#!/bin/bash

echo "*** somatic SNV by varScan ***"

TR=""
optTR=""
iniFile="`dirname $0`/../parameter/init_human.hg38.sh"
source $iniFile
_refData=$refData
###_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/bwa_0.7.12/human_g1k_v37_decoy.fasta"
#_refData="/workspace/zhengliangtao/00.database/broad/cellranger/human/genome.fa"
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

##source $iniFile
refData=$_refData
REF=$_refData

##module load samtools/1.14
##module load varScan/2.4.2

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir

## outDir sampleID normalBam tumorBam optTR refData
run_varScan()
{
    local outDir=$1
    local sampleID=$2
    local normalBam=$3
    local tumorBam=$4
    local optTR=$5
    local refData=$6
    #optTR="-l /workspace/zhengliangtao/00.database/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.exon.chr.ext5.bed"
    #refData="/workspace/zhengliangtao/00.database/broad/cellranger/human/genome.fa"
    REF=$refData
    #tumorFreq=0.1
    ### for scRNA-seq, normalFreq=0.001
    #normalFreq=0.001
    mkdir -p $outDir

    samtools mpileup -B -q 1 $optTR -f $refData $normalBam $tumorBam \
        | java -Xmx8G -jar $VarScanJAR somatic --mpileup 1 \
            --min-coverage-normal 8 \
            --min-coverage-tumor 8 \
            --p-value 0.1 \
            --min-var-freq $tumorFreq \
            --output-snp $outDir/$sampleID.varScan.snp \
            --output-indel $outDir/$sampleID.varScan.indel \
            --output-vcf

    java -Xmx4G -jar $VarScanJAR processSomatic \
        $outDir/$sampleID.varScan.snp.vcf \
        --min-tumor-freq $tumorFreq \
        --max-normal-freq $normalFreq

    java -Xmx4G -jar $VarScanJAR processSomatic \
        $outDir/$sampleID.varScan.indel.vcf \
        --min-tumor-freq $tumorFreq \
        --max-normal-freq $normalFreq

}

### call ###
if [ ! -f "$outDir/$sampleID.varScan.snp.Somatic.hc.vcf" ]; then
    #samtools index $tumorBam
    #samtools index $normalBam
    run_varScan $outDir $sampleID $normalBam $tumorBam "$optTR" "$refData"
fi

### fpfilter
## invcf tumorBam refData
filterSomatiMut()
{
    local invcf=$1
    local tumorBam=$2
    local refData=$3
    local REF=$refData

#    local vType=$2
#    if [ "$vtype" = "snp" ]; then
#        _opt=""
#    else
#        _opt="--insertion-centric"
#    fi

    awk -v OFS="\t" '!/^#/ {print $1,$2,$2+length($4)-1}' $invcf > ${invcf/.vcf/.list}

    bam-readcount $tumorBam -w 3 -q 1 -b 13 -f $refData -l ${invcf/.vcf/.list} >  ${invcf/.vcf/.readcount}

    ### note: --max-mapqual-diff [ default 50] --min-var-dist3 [default 0.1]
    ### use "conventional" parameters
    java -Xmx4G -jar $VarScanJAR fpfilter \
        $invcf \
        ${invcf/.vcf/.readcount} \
        --output-file ${invcf/.vcf/.fpfilter.pass.vcf} \
        --filtered-file ${invcf/.vcf/.fpfilter.fail.vcf} \
        --min-var-readpos 0.1 \
        --min-var-freq 0.05 \
        --min-var-count 4 \
        --min-var-count-lc 4 \
        --min-strand-reads 4 \
        --min-strandedness 0.01 \
        --max-mmqs-diff 50 \
        --max-mapqual-diff 30 \
        --max-rl-diff 0.25 \
        --min-var-dist3 0.2 \
        --max-var-mmqs 100

    ### donot norm vcf before bam-readcount !!!!
    bgzip -f ${invcf/.vcf/.fpfilter.pass.vcf}
    tabix -f -p vcf ${invcf/.vcf/.fpfilter.pass.vcf.gz}

    bcftools filter ${invcf/.vcf/.fpfilter.pass.vcf.gz} \
        -i 'FILTER=="PASS"||FILTER=="."' \
        | bcftools norm -m-both \
        | bcftools norm -f $REF \
        > ${invcf/.vcf/.fpfilter.pass.norm.vcf}

}

filterSomatiMut $outDir/$sampleID.varScan.snp.Somatic.hc.vcf $tumorBam $refData
filterSomatiMut $outDir/$sampleID.varScan.indel.Somatic.hc.vcf $tumorBam $refData

### annovarFile  annVer HumanDB
run_annovar()
{
    local annovarFile=$1
    local annVer=$2
    local HumanDB=$3
    #local annVer="hg38"
    #HumanDB="/workspace/zhengliangtao/01.bin/ANNOVAR/annovar/humandb"

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

run_annovar $outDir/$sampleID.varScan.snp.Somatic.hc.fpfilter.pass.norm.vcf "hg38" $HumanDB
run_annovar $outDir/$sampleID.varScan.indel.Somatic.hc.fpfilter.pass.norm.vcf "hg38" $HumanDB


echo end at: `date`
echo "*** Finished somatic SNV by varScan ***"

######################### old below ############################
######awk -v OFS="\t" '!/^#/ {print $1,$2,$2+length($4)-1}' $outDir/$sampleID.varScan.snp.Somatic.hc.vcf > $outDir/$sampleID.varScan.snp.Somatic.hc.list
######bam-readcount $tumorBam -w 3 -q 1 -b 13 -f $refData -l $outDir/$sampleID.varScan.snp.Somatic.hc.list > $outDir/$sampleID.varScan.snp.Somatic.hc.readcount
######$varScanDIR/fpfilter_vcf.pl $outDir/$sampleID.varScan.snp.Somatic.hc.vcf $outDir/$sampleID.varScan.snp.Somatic.hc.readcount --output-basename $outDir/$sampleID.varScan.snp.Somatic.hc.fpfilter
######awk '/^#/' $outDir/$sampleID.varScan.snp.Somatic.hc.vcf > $outDir/$sampleID.varScan.call.Somatic.vcf
######cat $outDir/$sampleID.varScan.snp.Somatic.hc.fpfilter.pass >> $outDir/$sampleID.varScan.call.Somatic.vcf
######
######## germline
######awk '!/^#/ {print $1,$2,$2}' $outDir/$sampleID.varScan.snp.Germline.hc.vcf > $outDir/$sampleID.varScan.snp.Germline.hc.list
######bam-readcount $tumorBam -w 3 -q 1 -b 13 -f $refData -l $outDir/$sampleID.varScan.snp.Germline.hc.list > $outDir/$sampleID.varScan.snp.Germline.hc.readcount
######$varScanDIR/fpfilter_vcf.pl $outDir/$sampleID.varScan.snp.Germline.hc.vcf $outDir/$sampleID.varScan.snp.Germline.hc.readcount --output-basename $outDir/$sampleID.varScan.snp.Germline.hc.fpfilter
######awk '/^#/' $outDir/$sampleID.varScan.snp.Germline.hc.vcf > $outDir/$sampleID.varScan.call.Germline.vcf
######cat $outDir/$sampleID.varScan.snp.Germline.hc.fpfilter.pass >> $outDir/$sampleID.varScan.call.Germline.vcf
######
######## LOH
######awk '!/^#/ {print $1,$2,$2}' $outDir/$sampleID.varScan.snp.LOH.hc.vcf > $outDir/$sampleID.varScan.snp.LOH.hc.list
######bam-readcount $tumorBam -w 3 -q 1 -b 13 -f $refData -l $outDir/$sampleID.varScan.snp.LOH.hc.list > $outDir/$sampleID.varScan.snp.LOH.hc.readcount
######$varScanDIR/fpfilter_vcf.pl $outDir/$sampleID.varScan.snp.LOH.hc.vcf $outDir/$sampleID.varScan.snp.LOH.hc.readcount --output-basename $outDir/$sampleID.varScan.snp.LOH.hc.fpfilter
######awk '/^#/' $outDir/$sampleID.varScan.snp.LOH.hc.vcf > $outDir/$sampleID.varScan.call.LOH.vcf
######cat $outDir/$sampleID.varScan.snp.LOH.hc.fpfilter.pass >> $outDir/$sampleID.varScan.call.LOH.vcf

