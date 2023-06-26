#!/bin/bash -eu


sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optP=8
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
	p)
		optP=$OPTARG
		;;
	m)
		optMode=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-p num. of CPUs. default 8] [-m mode. one of 'STAR', and ''. default ''] <sampleID> <fq1/bam> <fq2/-> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
    echo "Usage: $0 [-c iniFile] [-p num. of CPUs. default 8] [-m mode. one of 'STAR', and ''. default ''] <sampleID> <fq1/bam> <fq2/-> <outDir>"
	exit 1
fi

#source $iniFile

export KRAKEN2_DB_PATH="/lustre1/zeminz_pkuhpc/00.database/kraken/kraken2"
STAR_GenomeDir="/lustre1/zeminz_pkuhpc/00.database/gencode/release_23/vn/star"
gene_model_file="/lustre1/zeminz_pkuhpc/00.database/gencode/release_23/vn/gencode.v23.annotation.gtf"

export MODULESHOME=/usr/share/Modules
. /usr/share/Modules/init/bash
export MODULEPATH="/lustre1/zeminz_pkuhpc/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
module load blastPlus/2.2.30
module load kraken2/unknown
###module load samtools/0.1.19
module load samtools/1.14
module load picard/1.130
module load STAR/2.7.9a

sampleID=$1
fq1=$2
fq2=$3
outDir=$4

### for most analysis, parallel mode is not supported currently 
#optNT=1
#
#optL=''
#optM=`echo "scale=0;$optM/1.5" | bc`
#max_reads=`echo 250000*$optM | bc`

echo begin at: `date` 

mkdir -p $outDir

newFQ1=$fq1
newFQ2=$fq2
data_format="fastq"
bam_file=""

if [[ $fq1 =~ \.bam$ ]] && [[ $fq2 == "-" ]];then
    echo "input a bam file"
    bam_file=$fq1
    data_format="bam"
else
    echo "input two fastq files"
fi

if [[ $data_format == "bam" ]];then

    samtools view -f 12 --bam --output  $outDir/$sampleID.unmapped.PE.bam $bam_file && \
    samtools sort -n -m 4G -o $outDir/$sampleID.unmapped.PE.sortName.bam $outDir/$sampleID.unmapped.PE.bam && \
    rm $outDir/$sampleID.unmapped.PE.bam && \
    samtools fastq -O -1 $outDir/$sampleID.unmapped.1.fq.gz -2 $outDir/$sampleID.unmapped.2.fq.gz $outDir/$sampleID.unmapped.PE.sortName.bam && \
    rm $outDir/$sampleID.unmapped.PE.sortName.bam
    
    newFQ1=$outDir/$sampleID.unmapped.1.fq.gz
    newFQ2=$outDir/$sampleID.unmapped.2.fq.gz

fi

if [[ $optMode == "STAR" ]];then
    if [ ! -f "$outDir/$sampleID.Aligned.out.bam" ];then
        STAR \
            --runThreadN $optP \
            --genomeDir $STAR_GenomeDir \
            --readFilesIn $fq1 $fq2 \
            --readFilesCommand zcat \
            --outFileNamePrefix $outDir/$sampleID. \
            --quantMode GeneCounts \
            --twopassMode Basic \
            --outSAMattrRGline ID:$sampleID CN:BIOPIC LB:$sampleID PL:illumina PU:$sampleID SM:$sampleID \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within KeepPairs \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --chimSegmentReadGapMax 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outSAMstrandField intronMotif \
            --chimOutJunctionFormat 1  
    fi

    samtools view -f 12 --bam --output  $outDir/$sampleID.unmapped.PE.bam $outDir/$sampleID.Aligned.out.bam && \
    samtools sort -n -m 4G -o $outDir/$sampleID.unmapped.PE.sortName.bam $outDir/$sampleID.unmapped.PE.bam && \
    rm $outDir/$sampleID.unmapped.PE.bam && \
    samtools fastq -O -1 $outDir/$sampleID.unmapped.1.fq.gz -2 $outDir/$sampleID.unmapped.2.fq.gz $outDir/$sampleID.unmapped.PE.sortName.bam && \
    rm $outDir/$sampleID.unmapped.PE.sortName.bam && \
    rm -r $outDir/$sampleID._STAR* && \
    rm $outDir/$sampleID.Aligned.out.bam

    newFQ1=$outDir/$sampleID.unmapped.1.fq.gz
    newFQ2=$outDir/$sampleID.unmapped.2.fq.gz

fi

kraken2 --db STDL --threads $optP --gzip-compressed --paired \
    --report-zero-counts \
    --use-names \
    --use-mpa-style \
    --minimum-base-quality 1 \
    --report $outDir/$sampleID.kraken2.report.out \
    --output $outDir/$sampleID.kraken2.detail.out \
    $newFQ1 $newFQ2

echo end at: `date` 
