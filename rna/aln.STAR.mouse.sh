#!/bin/bash -eu


sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optM=8
optMode=""

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
		optMode=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-m mode, default ''] <sampleID> <fq1> <fq2> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-m mode, default ''] <sampleID> <fq1> <fq2> <outDir>"
	exit 1
fi

source $iniFile

STAR_GenomeDir=/WPSnew/zhenglt/work/proj_fh/data/refdata-cellranger-mm10-1.2.0/star
gene_model_file=/WPSnew/zhenglt/work/proj_fh/data/refdata-cellranger-mm10-1.2.0/genes/genes.gtf

export PICARD="/WPSnew/zhenglt/01.bin/var/picard/picard.2.18.9"
module load STAR/2.6.0c
module load samtools/1.8
###module unload java/1.7.0_79
module load java/1.8.0_171
module load subread/1.6.2
###module unload gatk/3.3-0
###module load gatk/3.8-0
###module load trimmomatic/0.33

export REF=/WPSnew/zhenglt/work/proj_fh/data/refdata-cellranger-mm10-1.2.0/fasta/genome.fa
### with ERCC
### without ERCC
#export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/hg19.order.after.gsnap.fa
####knownSites=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/dbsnp_138.hg19.reorder.vcf

sampleID=$1
fq1=$2
fq2=$3
outDir=$4
#REF=$5

optL=''
optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`

mkdir -p $outDir

newFQ1=$fq1
if [ "$fq2" == "-" ];then
    newFQ2=""
else
    newFQ2=$fq2
fi
echo newFQ1: $newFQ1
echo newFQ2: $newFQ2
len_trim_to=101

#if [ "$optMode" == "trim" ];then
#	java -Xmx${optM}g -jar $TRIMMOMATIC_DIR/trimmomatic-0.33.jar PE \
#		-threads 8 \
#		$fq1 \
#		$fq2 \
#		$outDir/$sampleID.clean.P.R1.fq.gz \
#		$outDir/$sampleID.clean.UP.R1.fq.gz \
#		$outDir/$sampleID.clean.P.R2.fq.gz \
#		$outDir/$sampleID.clean.UP.R2.fq.gz \
#		ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/NexteraPE-PE.fa:2:30:10 \
#		CROP:$len_trim_to \
#		TRAILING:3 \
#		MAXINFO:50:0.25 \
#		MINLEN:36 \
#		TOPHRED33 
#	echo "trimmomatic done. (" `date` ")"
#    newFQ1=$outDir/$sampleID.clean.P.R1.fq.gz
#    newFQ2=$outDir/$sampleID.clean.P.R2.fq.gz
#fi

STAR \
    --runThreadN 8 \
    --genomeDir $STAR_GenomeDir \
    --readFilesIn $newFQ1 $newFQ2 \
    --readFilesCommand zcat \
    --outFileNamePrefix $outDir/$sampleID. \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --outSAMattrRGline ID:$sampleID CN:BIOPIC LB:$sampleID PL:illumina PU:$sampleID SM:$sampleID \
    --outSAMtype BAM SortedByCoordinate
samtools index $outDir/${sampleID}.Aligned.sortedByCoord.out.bam

echo "... using $max_reads reads in memory (parameter optM: $optM*1.5)"
java -Xmx${optM}g -jar $PICARD/picard.jar SortSam \
		I=$outDir/${sampleID}.Aligned.sortedByCoord.out.bam \
		O=$outDir/${sampleID}.Aligned.sortedByRName.out.bam \
		MAX_RECORDS_IN_RAM=$max_reads \
		TMP_DIR=$outDir \
		SO=queryname \
		VALIDATION_STRINGENCY=SILENT
#samtools index $outDir/${sampleID}.Aligned.sortedByRName.out.bam

featureCounts -p -T 2 \
    -a $gene_model_file \
    -o $outDir/$sampleID.subread.exp \
    $outDir/${sampleID}.Aligned.sortedByRName.out.bam
    ##$outDir/${sampleID}.Aligned.sortedByCoord.out.bam

featureCounts -p -T 2 -O \
    -a $gene_model_file \
    -o $outDir/$sampleID.subread.exp.optO \
    $outDir/${sampleID}.Aligned.sortedByRName.out.bam
    ##$outDir/${sampleID}.Aligned.sortedByCoord.out.bam

rm $outDir/${sampleID}.Aligned.sortedByRName.out.bam    

echo end at: `date` 
