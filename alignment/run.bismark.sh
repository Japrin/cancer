#!/bin/bash


iniFile="/WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/parameter/init_human.sh"
_refData="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optE=""
optM=16

while getopts c:m:e: opt
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
	e)
		optE=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-e etra parameters] [-m memory, default 16(GB)] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-c iniFile] [-e etra parameters] [-m memory, default 16(GB)] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi

echo begin at: `date`

source $iniFile
refData=$_refData
REF=$_refData

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
lib=$5
lane=$6

mkdir -p $outDir

genome_folder=`dirname $REF`/bismark

mySortBam()
{
	inbam=$1
	outbam=$2
	SO=$3
	#SO="coordinate"
	optM=`echo "scale=0;$optM/1.5" | bc`
	max_reads=`echo 250000*$optM | bc`
	echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
	java -Xmx${optM}g -jar $picardDIR/SortSam.jar \
		TMP_DIR=$outDir \
		I=$inbam \
		O=$outbam \
		MAX_RECORDS_IN_RAM=$max_reads \
		TMP_DIR=$outDir \
		SO=$SO \
		VALIDATION_STRINGENCY=SILENT
	samtools index $outbam
}

myDeDup()
{
	inbam=$1
	outbam=$2
	rmdup="true"
	optM=`echo "scale=0;$optM/1.5" | bc`
	max_reads=`echo 250000*$optM | bc`
	echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
	java -Xmx${optM}g -jar $PICARD/MarkDuplicates.jar \
		TMP_DIR=$outDir \
		I=$inbam \
		O=$outbam \
		M=${outbam/.bam/.metrics} \
		VALIDATION_STRINGENCY=SILENT \
		ASSUME_SORTED=true \
		REMOVE_DUPLICATES=$rmdup \
		MAX_RECORDS_IN_RAM=$max_reads
	samtools index $outbam
}

if [ "$fq2" != '-' ]; then
	bismark $optE --unmapped --ambiguous --output_dir $outDir --bowtie2 -p 4 --temp_dir $outDir --non_bs_mm --gzip --bam --basename $sampleID.$lib $genome_folder -1 $fq1 -2 $fq2
	mySortBam $outDir/$sampleID.${lib}_pe.bam $outDir/$sampleID.${lib}_pe.sort.bam "coordinate"
	myDeDup   $outDir/$sampleID.${lib}_pe.sort.bam $outDir/$sampleID.${lib}_pe.sort.rmdup.bam
	mySortBam $outDir/$sampleID.${lib}_pe.sort.rmdup.bam $outDir/$sampleID.${lib}_pe.sort.rmdup.resort.bam "queryname"
	finalBam=$outDir/$sampleID.${lib}_pe.sort.rmdup.resort.bam
	bismark_methylation_extractor \
		-p --comprehensive --report --output $outDir --gzip \
		--bedGraph --zero_based --buffer_size 8G \
		--cytosine_report --genome_folder $genome_folder \
		--no_overlap  $finalBam
else
	bismark $optE --unmapped --ambiguous --output_dir $outDir --bowtie2 -p 4 --temp_dir $outDir --non_bs_mm --gzip --bam --basename $sampleID.$lib $genome_folder $fq1
	mySortBam $outDir/$sampleID.${lib}_pe.bam $outDir/$sampleID.${lib}_pe.sort.bam "coordinate"
	myDeDup   $outDir/$sampleID.${lib}_pe.sort.bam $outDir/$sampleID.${lib}_pe.sort.rmdup.bam
	finalBam=$outDir/$sampleID.${lib}_pe.sort.rmdup.bam
	bismark_methylation_extractor \
		-s --comprehensive --report --output $outDir --gzip \
		--bedGraph --zero_based --buffer_size 8G \
		--cytosine_report --genome_folder $genome_folder \
		$finalBam
fi



echo end at: `date`
