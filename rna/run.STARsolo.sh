#!/bin/bash -eu


sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optM=30
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

#source $iniFile

STAR_GenomeDir=/workspace/zhengliangtao/00.database/genome/STARsolo/human/star
gene_model_file=/workspace/zhengliangtao/00.database/genome/STARsolo/human/genes.gtf
#STAR_FUSION_DB_Dir=/WPSnew/zhenglt/00.database/tools/star-fusion/GRCh37_v24_mybuild/ctat_genome_lib_build_dir

module load STAR/2.7.9a
#module load java/1.8.0_171
#module load subread/1.6.3
#module load trimmomatic/0.33
#module load blast/2.8.1+
#module load STAR-Fusion/1.5.0

### with ERCC
#####export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/withERCC.v1/hg19.order.after.gsnap.fa
### without ERCC
#####export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/hg19.order.after.gsnap.fa
#####knownSites=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/dbsnp_138.hg19.reorder.vcf

sampleID=$1
fq1=$2
fq2=$3
outDir=$4
####REF=$5

### for most analysis, parallel mode is not supported currently 
optNT=1
###optNT=8

optL=''
optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`

mkdir -p $outDir
cd $outDir

newFQ1=$fq1
newFQ2=$fq2
len_trim_to=101

### Importantly, in the --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read, i.e.
### --readFilesIn cDNAfragmentSequence.fastq.gz CellBarcodeUMIsequence.fastq.gz
###
### Matching CellRanger 3.x.x results
###    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
###
### Matching CellRanger 4.x.x and 5.x.x results
###    --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
###    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
###
### --clip5pNbases and --soloBarcodeMate 1 cannot be used with --clipAdapterType CellRanger4. The latter option is specific to 10X adapter sequence and should only be used for 10X data.

echo begin at: `date`

    ##--soloBarcodeMate 1   --clip5pNbases 39 0 \

STAR \
    --runThreadN 8 \
    --genomeDir $STAR_GenomeDir \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple   --soloCBstart 1   --soloCBlen 16   --soloUMIstart 17   --soloUMIlen 10 \
    --soloBarcodeReadLength 0 \
    --readFilesIn $newFQ2 $newFQ1 \
    --soloCBwhitelist $sDir/737K-august-2016.txt \
    --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
    --soloCellFilter  EmptyDrops_CR \
    --soloFeatures Gene Velocyto \
    --soloMultiMappers EM \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $outDir/$sampleID. \
    --outSAMattrRGline ID:$sampleID CN:BIOPIC LB:$sampleID PL:illumina PU:$sampleID SM:$sampleID \
    --outReadsUnmapped Fastx \
    --outSAMunmapped Within \
    --twopassMode Basic \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --alignSJDBoverhangMin 10 \
    --alignMatesGapMax 100000 \
    --alignIntronMax 100000 \
    --chimSegmentReadGapMax 3 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --outSAMstrandField intronMotif \
    --chimOutJunctionFormat 1
    #--quantMode GeneCounts \

####### gene expression

echo end at: `date` 
