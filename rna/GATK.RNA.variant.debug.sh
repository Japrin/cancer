#!/bin/bash -eu

echo "*** Realigning targeted regions ***"

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
optM=8
optNT=1
optNCT=4

while getopts c:m: opt
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
	m)
		optM=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] <sampleID> <inbam> <outDir> <REF>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] <sampleID> <inbam> <outDir> <REF>"
	exit 1
fi

source $iniFile

module unload java/1.7.0_79
module load java/1.8.0_144
module unload gatk/3.3-0
module load gatk/3.8-0
### with ERCC
#export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/withERCC.v1/hg19.order.after.gsnap.fa
### without ERCC
#export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/hg19.order.after.gsnap.fa
knownSites=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/dbsnp_138.hg19.reorder.vcf

sampleID=$1
inBam=$2
outDir=$3
REF=$4

mkdir -p $outDir

optL=''

outBam=$outDir/$sampleID.GATK.RNA.ready.bam
optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`

if [ -f $outBam ]
then
    echo "## $outBam exists, skip this step"
	exit 0
fi

if [ ! -f $inBam.bai ]
then
	samtools index $inBam
fi

#echo ">>> Marking duplicates"
#this_bam=$outDir/$sampleID.rmDup.bam
#if [ ! -f $this_bam ]
#then
#    echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
#    java -Xmx${optM}g -jar $PICARD/picard.jar MarkDuplicates \
#    	TMP_DIR=$outDir \
#    	I=$inBam \
#    	O=$this_bam \
#    	M=${this_bam%.bam}.metrics \
#    	VALIDATION_STRINGENCY=SILENT \
#    	ASSUME_SORTED=true \
#    	REMOVE_DUPLICATES=false \
#    	MAX_RECORDS_IN_RAM=$max_reads
#    
#    rm -f ${this_bam%.bam}.bai
#    samtools index $this_bam
#fi
#
#echo ">>> splitNTrim"
#optM=`echo "scale=0;$optM/1.5" | bc`
#inBam=$this_bam
#this_bam=$outDir/$sampleID.GATK.RNA.splitNCigar.bam
#if [ ! -f $this_bam.bai ]
#then
##    java -Xmx${optM}g -jar $PICARD/picard.jar FixMateInformation \
##    	TMP_DIR=$outDir \
##    	I=$inBam \
##    	O=$this_bam.fix.bam \
##    	VALIDATION_STRINGENCY=SILENT \
##    	ASSUME_SORTED=true \
##    	MAX_RECORDS_IN_RAM=$max_reads
##    samtools index $this_bam.fix.bam
#
#    java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
#        -T SplitNCigarReads \
#        -R $REF \
#        -I $inBam \
#        -o $this_bam \
#        -U ALLOW_N_CIGAR_READS \
#	    -nt $optNT 
#        ##-fixNDN
#	#-et NO_ET \
#	#-K $GATKKey
#    ## for gsnap result, not nessesary:
#    ##-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
#    samtools index $this_bam
#fi
#
#rm $inBam
#rm $inBam.bai
#
#echo ">>> Determining (small) suspicious intervals which are likely in need of realignment"
#java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
#	-T RealignerTargetCreator \
#	-I $outDir/$sampleID.GATK.RNA.splitNCigar.bam \
#	-R $REF \
#	-o $outDir/$sampleID.GATK.RNA.realn.intervals $optL \
#	-nt $optNT \
#	-rf BadCigar
#echo ">>> Running the realigner over the targeted intervals"
#java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
#	-T IndelRealigner \
#	-I $outDir/$sampleID.GATK.RNA.splitNCigar.bam \
#	-R $REF \
#	-o $outDir/$sampleID.GATK.RNA.realn.bam \
#	-targetIntervals $outDir/$sampleID.GATK.RNA.realn.intervals \
#	-LOD 5 $optL \
#	-nt $optNT \
#	-rf BadCigar
#samtools index $outDir/$sampleID.GATK.RNA.realn.bam
#
#java -Xmx${optM}g -jar $PICARD/picard.jar AddOrReplaceReadGroups \
#    TMP_DIR=$outDir \
#    I=$outDir/$sampleID.GATK.RNA.realn.bam \
#    O=$outDir/$sampleID.GATK.RNA.realn.RG.bam \
#    VALIDATION_STRINGENCY=SILENT \
#    MAX_RECORDS_IN_RAM=$max_reads \
#    RGID=$sampleID RGLB=$sampleID RGPL=illumina RGPU=$sampleID RGSM=$sampleID
#samtools index $outDir/$sampleID.GATK.RNA.realn.RG.bam
#rm $outDir/$sampleID.GATK.RNA.realn.bam
#
#echo ">>> Counting covariates"
#optM=`echo "scale=0;$optM/1.2" | bc`
#java -Djava.io.tmpdir=$outDir -Xms${optM}g -Xmx${optM}g -jar $GATK/GenomeAnalysisTK.jar \
#	-T BaseRecalibrator \
#	-R $REF \
#	--knownSites $knownSites \
#	--disable_indel_quals \
#	-I $outDir/$sampleID.GATK.RNA.realn.RG.bam \
#	-o $outDir/$sampleID.GATK.RNA.recal.grp \
#	-cov ReadGroupCovariate \
#	-cov QualityScoreCovariate \
#	-cov CycleCovariate \
#	-cov ContextCovariate \
#	-rf BadCigar \
#    --validation_strictness SILENT
#
#echo ">> Table recalibration"
#if [ "`grep -v '#' ${o/.bam/.grp} | grep -v "EOF" | wc -l`" = "1" ]
#then
#    ln $outDir/$sampleID.GATK.RNA.realn.RG.bam $outDir/$sampleID.GATK.RNA.recal.bam
#else
#	java -Djava.io.tmpdir=$outDir -Xms${optM}g -Xmx${optM}g -jar $GATK/GenomeAnalysisTK.jar \
#	-T PrintReads \
#	-R $REF \
#	-I $outDir/$sampleID.GATK.RNA.realn.RG.bam \
#	-o $outDir/$sampleID.GATK.RNA.recal.bam \
#	-BQSR $outDir/$sampleID.GATK.RNA.recal.grp \
#	--disable_indel_quals \
#	-EOQ \
#	-rf BadCigar
#fi
#samtools index $outDir/$sampleID.GATK.RNA.recal.bam

### remove ERCC
samtools view -h $outDir/$sampleID.GATK.RNA.recal.bam | awk -F"\t" '!/SN:ERCC/ && $3!~/^ERCC/ && $7!~/^ERCC/' \
    | samtools view -b -o $outDir/$sampleID.GATK.RNA.recal.fltERCC.bam
samtools index $outDir/$sampleID.GATK.RNA.recal.fltERCC.bam

java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $REF \
    -I $outDir/$sampleID.GATK.RNA.recal.bam \
    -dontUseSoftClippedBases \
    -stand_call_conf 20.0 \
	-nt $optNT \
    --num_cpu_threads_per_data_thread 4 \
    -o $outDir/$sampleID.GATK.RNA.recal.var.vcf

java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $REF \
    -V $outDir/$sampleID.GATK.RNA.recal.var.vcf \
    -window 35 -cluster 3 \
    -filterName FS -filter "FS > 30.0" \
    -filterName QD -filter "QD < 2.0" \
    -o $outDir/$sampleID.GATK.RNA.recal.var.flt.vcf


echo end at: `date` 
