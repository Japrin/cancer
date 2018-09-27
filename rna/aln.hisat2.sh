#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optT=8
optA=""
optB=""
optM=16
optS="human"

while getopts c:t:f:a:b:m:s: opt
do
	case $opt in 
	t)
		optT=$OPTARG
		;;
    a)
        #optA="--phred64"
        optA="$OPTARG"
        ;;
    s)
        optS=$OPTARG
        ;;
	b)
		optB="--upto $OPTARG"
		;;
	m)
		optM=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-t threads, default 8] [-s species, default human] [-a quality mode, \"--phred33\" or \"--phred64\" (default \"--phred33\")] [-b upto] [-m memrory(GB), default 16] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-t threads, default 8] [-s species, default human] [-a quality mode, \"--phred33\" or \"--phred64\" (default \"--phred33\")] [-b upto] [-m memrory(GB), default 16] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi

source $iniFile
module load hisat2/2.0.5
module load bamUtil/1.0.14
module load trimmomatic/0.33
module load RNA-SeQC/1.1.8
####module load bwa/0.7.12
##module load bwa/0.7.17
#module load picard/1.130
module load picard/2.18.9
### pipeline need samtools 1.2; samtools 1.3.1 has some problem
module load samtools/1.2

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
lib=$5
lane=$6
TMPDIR="/WPSnew/zhenglt/tmp"

echo begin at: `date`
echo "*** Aligning reads ***"

mkdir -p $outDir
#export HISAT_INDEXES="/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/hisat"
#INDEX="$HISAT_INDEXES/hg19"
#KNOWN_SPLICE="$HISAT_INDEXES/hg19_knownGene.splicesites.txt"
#RRNA_INDEX="$HISAT_INDEXES/Rfam.human.rRNA"
if [ "$optS" == "human" ];then
    export HISAT_INDEXES="/DBS/DB_temp/zhangLab/gencode/release_19/hisat"
    INDEX="$HISAT_INDEXES/GRCh37.p13.genome"
    KNOWN_SPLICE="$HISAT_INDEXES/gencode.v19.splicesites.txt"
    RRNA_INDEX="$HISAT_INDEXES/human_all_rRNA"
    GENE_MODEL="$HISAT_INDEXES/gencode.v19.annotation.gtf"
    #GCFile="$HISAT_INDEXES/gencode.v19.annotation.GC"
    REF="$INDEX.fa"
elif [ "$optS" == "mouse" ];then
    export HISAT_INDEXES="/WPSnew/zhenglt/00.database/ensemble/mybuild/hisat/mouse/rel89/hisat2"
    INDEX="$HISAT_INDEXES/GRCm38.rel89"
    KNOWN_SPLICE="$HISAT_INDEXES/GRCm38.rel89.splicesites.txt"
    GENE_MODEL="$HISAT_INDEXES/GRCm38.rel89.gtf"
    REF="$INDEXES.fa"
    RRNA_INDEX="$HISAT_INDEXES/mouse_all_rRNA"
fi

rm -f $outDir/$sampleID.noRRNA.R1.fq.gz $outDir/$sampleID.noRRNA.R2.fq.gz 
#mkfifo $outDir/$sampleID.noRRNA.R1.fq
#mkfifo $outDir/$sampleID.noRRNA.R2.fq
rm -f $outDir/$sampleID.clean.P.R1.fq.gz $outDir/$sampleID.clean.P.R2.fq.gz
#mkfifo $outDir/$sampleID.clean.P.R1.fq
#mkfifo $outDir/$sampleID.clean.P.R2.fq
#echo $outDir/$sampleID.clean.UP.R1.fq
#echo $outDir/$sampleID.clean.UP.R2.fq


if [ ! -f "$outDir/$sampleID.hisat.hit.sort.bam" ] && [ ! -f "$outDir/$sampleID.hisat.hit.bam" ];then
	### filter low quality reads
		##ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-PE-2.fa:2:30:10 \
	java -Xmx${optM}g -jar $TRIMMOMATIC_DIR/trimmomatic-0.33.jar PE \
		-threads 8 \
		$fq1 \
		$fq2 \
		$outDir/$sampleID.clean.P.R1.fq.gz \
		$outDir/$sampleID.clean.UP.R1.fq.gz \
		$outDir/$sampleID.clean.P.R2.fq.gz \
		$outDir/$sampleID.clean.UP.R2.fq.gz \
		ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/NexteraPE-PE.fa:2:30:10 \
		CROP:75 \
		TRAILING:3 \
		MAXINFO:50:0.25 \
		MINLEN:36 \
		TOPHRED33 
	echo "trimmomatic done. (" `date` ")"
	
	### filter contaminated rRNA
	hisat2  -x $RRNA_INDEX -1 $outDir/$sampleID.clean.P.R1.fq.gz -2 $outDir/$sampleID.clean.P.R2.fq.gz \
		--minins 0 \
		--maxins 1000 \
		--fr \
		-t \
		--met-file "$outDir/$sampleID.hisat.rrna.met.txt" \
		--met 10 \
		--rg-id $sampleID \
		--rg "SM:$sampleID" \
		--rg "LB:$lib" \
		--rg "PU:$lane" \
		--threads $optT $optA $optB \
		| samtools view -f 0xc - \
		| bam bam2FastQ --readName --in - --firstOut $outDir/$sampleID.noRRNA.R1.fq.gz --secondOut $outDir/$sampleID.noRRNA.R2.fq.gz
	
	rm $outDir/$sampleID.clean.P.R1.fq.gz $outDir/$sampleID.clean.P.R2.fq.gz
	echo "filtering rRNA done. (" `date` ")"
	### hisat alignment
	hisat2  -x $INDEX -1 $outDir/$sampleID.noRRNA.R1.fq.gz -2 $outDir/$sampleID.noRRNA.R2.fq.gz \
		--known-splicesite-infile $KNOWN_SPLICE \
		--novel-splicesite-outfile "$outDir/$sampleID.hisat.novel-splicesite.txt" \
		--minins 0 \
		--maxins 1000 \
		--fr \
		-t \
		--met-file "$outDir/$sampleID.hisat.met.txt" \
		--met 10 \
		--rg-id $sampleID \
		--rg "SM:$sampleID" \
		--rg "LB:$lib" \
		--rg "PU:$lane" \
		--threads $optT $optA \
		| samtools view -bS - > $outDir/$sampleID.hisat.hit.bam
	#rm $outDir/$sampleID.noRRNA.R1.fq.gz $outDir/$sampleID.noRRNA.R2.fq.gz
	echo "hisat mapping done. (" `date` ")"
	#wait
fi

if [ ! -f "$outDir/$sampleID.hisat.hit.sort.bam" ];then
	JM=`echo "scale=0;$optM/1.5" | bc`
	max_reads=`echo 250000*$JM | bc`
	echo "using $max_reads reads in memory (parameter optM: $optM)"
	java -Xmx${JM}g -jar $picardDIR/picard.jar SortSam \
			I=$outDir/$sampleID.hisat.hit.bam \
			O=$outDir/$sampleID.hisat.hit.sort.bam \
			MAX_RECORDS_IN_RAM=$max_reads \
			TMP_DIR=$outDir \
			SO=coordinate \
			VALIDATION_STRINGENCY=SILENT
	if [ -f "$outDir/$sampleID.hisat.hit.sort.bam" ];then
		#echo rm $outDir/$sampleID.hisat.hit.bam 
		rm $outDir/$sampleID.hisat.hit.bam 
	fi
	samtools index $outDir/$sampleID.hisat.hit.sort.bam
	echo "sorting bam done. (" `date` ")"
fi

if [ -f "$outDir/$sampleID.hisat.hit.sort.bam" ];then
	JM=`echo "scale=0;$optM/1.5" | bc`
    echo "java -Xmx${JM}g -jar $RNASeQC ..."
	#java -Xmx${JM}g -jar $RNASeQC \
	#	-n 1000 \
	#	-s "$sampleID|$outDir/$sampleID.hisat.hit.sort.bam|NA" \
	#	-t $GENE_MODEL \
	#	-r $REF \
	#	-o $outDir/RNA-SeQC \
	#	-strat gc -gc $GCFile \
	#	-BWArRNA "$RRNA_INDEX.fa"
	echo "RNASeQC done. (" `date` ")"
fi


echo end at: `date`
echo "*** Finished aligning reads ***"
