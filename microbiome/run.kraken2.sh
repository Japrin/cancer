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

export KRAKEN2_DB_PATH="/workspace/zhengliangtao/00.database/kraken2"
STAR_GenomeDir="/workspace/zhengliangtao/00.database/genome/STAR/build"
gene_model_file="/workspace/zhengliangtao/00.database/genome/STAR/genes.gtf"

#export MODULESHOME=/usr/share/Modules
#. /usr/share/Modules/init/bash
#export MODULEPATH="/lustre1/zeminz_pkuhpc/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
module load ncbi-blast/2.14.0
module load kraken2/github
###module load samtools/0.1.19
module load samtools/1.14
module load java/18.0.1.1
module load picard/3.0.0 
module load STAR/2.7.9a
module load fastp/0.23.4

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
cd $outDir

#newFQ1=$fq1
#newFQ2=$fq2
data_format="fastq"
bam_file=""

if [[ $fq1 =~ \.bam$ ]] && [[ $fq2 == "-" ]];then
    echo "input a bam file"
    bam_file=$fq1
    data_format="bam"
else
    echo "input two fastq files"
fi

newFQ1=$outDir/$sampleID.unmapped.1.fq.gz
newFQ2=$outDir/$sampleID.unmapped.2.fq.gz

if [[ $data_format == "bam" ]];then

    samtools view -f 12 --bam --output  $outDir/$sampleID.unmapped.PE.bam $bam_file && \
    samtools sort -n -m 4G -o $outDir/$sampleID.unmapped.PE.sortName.bam $outDir/$sampleID.unmapped.PE.bam && \
    rm $outDir/$sampleID.unmapped.PE.bam && \
    samtools fastq -O -1 $outDir/$sampleID.unmapped.00.R1.fq.gz -2 $outDir/$sampleID.unmapped.00.R2.fq.gz $outDir/$sampleID.unmapped.PE.sortName.bam && \
    rm $outDir/$sampleID.unmapped.PE.sortName.bam
    
    fastp -i $outDir/$sampleID.unmapped.00.R1.fq.gz -I $outDir/$sampleID.unmapped.00.R2.fq.gz \
          -o $newFQ1 -O $newFQ2
    
else
    fastp -i $fq1 -I $fq2 -o $newFQ1 -O $newFQ2
fi



if [[ $optMode == "STAR" ]];then
    if [ ! -f "$outDir/$sampleID.unmapped.00.R1.fq.gz" ];then
        if [ ! -f "$outDir/$sampleID.Aligned.out.bam" ];then
            STAR \
                --runThreadN $optP \
                --genomeDir $STAR_GenomeDir \
                --readFilesIn $newFQ1 $newFQ2 \
                --readFilesCommand zcat \
                --outFileNamePrefix $outDir/$sampleID. \
                --quantMode GeneCounts \
                --twopassMode None \
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
        samtools fastq -O -1 $outDir/$sampleID.unmapped.00.R1.fq.gz -2 $outDir/$sampleID.unmapped.00.R2.fq.gz $outDir/$sampleID.unmapped.PE.sortName.bam && \
        rm $outDir/$sampleID.unmapped.PE.sortName.bam && \
        rm -r $outDir/$sampleID._STAR* && \
        rm $outDir/$sampleID.Aligned.out.bam
    
        newFQ1=$outDir/$sampleID.unmapped.00.R1.fq.gz
        newFQ1=$outDir/$sampleID.unmapped.00.R2.fq.gz

    fi
fi


kraken2 --db STDL --threads $optP --gzip-compressed --paired \
    --report-zero-counts \
    --use-names \
    --use-mpa-style \
    --minimum-base-quality 1 \
    --report $outDir/$sampleID.kraken2.report.00.out \
    --output $outDir/$sampleID.kraken2.detail.00.out \
    $newFQ1 $newFQ2

gzip -v $outDir/$sampleID.kraken2.report.00.out
gzip -v $outDir/$sampleID.kraken2.detail.00.out

echo end at: `date` 
