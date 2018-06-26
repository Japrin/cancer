#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optBam=false
while getopts b opt
do
	case $opt in 
	b)	
		optBam=true
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-b bamfile] <outDir> <sampleID> <fq1> [fq2]"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))


if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-b bamfile] <outDir> <sampleID> <fq1> [fq2]"
	exit 1
fi

outDir=$1
sampleID=$2
fq1=$3
fq2=$4

optT=4

mkdir -p $outDir

source $iniFile
. /usr/share/Modules/init/bash
export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
module load java/1.8.0_112
module load bamUtil/1.0.12

mixcrBin="/Share/BP/zhenglt/01.bin/repSeq/mixcr-2.0.2/mixcr"
optM=6
echo begin at: `date`
echo $optBam
if $optBam;then
    max_reads=`echo 250000*$optM | bc`
    echo java -Xmx${optM}g -jar $picardDIR/picard.jar SortSam I=$fq1 O=$outDir/$sampleID.sort.name.bam MAX_RECORDS_IN_RAM=$max_reads TMP_DIR=$outDir SO=queryname VALIDATION_STRINGENCY=SILENT
    java -Xmx${optM}g -jar $picardDIR/picard.jar SortSam I=$fq1 O=$outDir/$sampleID.sort.name.bam MAX_RECORDS_IN_RAM=$max_reads TMP_DIR=$outDir SO=queryname VALIDATION_STRINGENCY=SILENT
    ### for STAR alignment
    samtools view -F 0x100 $outDir/$sampleID.sort.name.bam | bam bam2FastQ --readName --in - --firstOut $outDir/$sampleID.forMiXCR.R1.fq.gz --secondOut $outDir/$sampleID.forMiXCR.R2.fq.gz
    ######bam bam2FastQ --readName --in $outDir/$sampleID.sort.name.bam --firstOut $outDir/$sampleID.forMiXCR.R1.fq.gz --secondOut $outDir/$sampleID.forMiXCR.R2.fq.gz
    fq1=$outDir/$sampleID.forMiXCR.R1.fq.gz
    fq2=$outDir/$sampleID.forMiXCR.R2.fq.gz
    echo rm $outDir/$sampleID.sort.name.bam
    rm $outDir/$sampleID.sort.name.bam
else
    echo "nothing to do"
fi
echo fq1: $fq1 
echo fq2: $fq2

$mixcrBin align -f -p rna-seq -OallowPartialAlignments=true --save-description --save-reads -t $optT \
    -r $outDir/$sampleID.MiXCR.log.align.txt \
    $fq1 $fq2 \
    $outDir/$sampleID.MiXCR.alignments.vdjca
$mixcrBin assemblePartial -f -p \
    -r $outDir/$sampleID.MiXCR.log.assemble.txt \
    $outDir/$sampleID.MiXCR.alignments.vdjca \
    $outDir/$sampleID.MiXCR.alignmentsRescued_1.vdjca
$mixcrBin assemblePartial -f -p \
    -r $outDir/$sampleID.MiXCR.log.assemble.txt \
    $outDir/$sampleID.MiXCR.alignmentsRescued_1.vdjca \
    $outDir/$sampleID.MiXCR.alignmentsRescued_2.vdjca
$mixcrBin assemble -f -t $optT -OaddReadsCountOnClustering=true -ObadQualityThreshold=15 \
    -r $outDir/$sampleID.MiXCR.log.assembleClones.txt \
    $outDir/$sampleID.MiXCR.alignmentsRescued_2.vdjca \
    $outDir/$sampleID.MiXCR.clones.clns
$mixcrBin exportClones -f $outDir/$sampleID.MiXCR.clones.clns $outDir/$sampleID.MiXCR.clones.txt
$mixcrBin exportClonesPretty $outDir/$sampleID.MiXCR.clones.clns $outDir/$sampleID.MiXCR.clonesPretty.txt
$mixcrBin exportClones -f --chains TCR $outDir/$sampleID.MiXCR.clones.clns $outDir/$sampleID.MiXCR.TCR.clones.txt
$mixcrBin exportClonesPretty --chains TCR $outDir/$sampleID.MiXCR.clones.clns $outDir/$sampleID.MiXCR.TCR.clonesPretty.txt

if $optBam;then
    echo rm $fq1 $fq2
    rm $fq1 $fq2
fi

/Share/BP/zhenglt/02.pipeline/cancer/repseq/TCRasm.MiXCR.slim.pl \
    -s $sampleID \
    $outDir/$sampleID.MiXCR.TCR.clones.txt \
    $outDir/$sampleID.MiXCR.TCR.clones.slim.txt

echo end at: `date`
