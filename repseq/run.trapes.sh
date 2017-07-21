#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optT=8

optBam=false
while getopts t:b opt
do
	case $opt in 
    t)
        optT=$OPTARG
        ;;
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

#mkdir -p $outDir/input/$sampleID
#sampleID=SAMPLEID
#inDir=/Share/BP/zhenglt/01.bin/RepSeq/VDJPuzzle/Example
#outDir=/Share/BP/zhenglt/01.bin/RepSeq/VDJPuzzle/test

###source $iniFile

export MODULESHOME=/usr/share/Modules
. /usr/share/Modules/init/bash
export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
module load novomedSeq/1.0
module load bowtie2/2.2.3
##module load samtools/0.1.19
module load rsem/1.3.0
##module load igblast/1.4.0
##module load java/1.7.0_79
##module load trinity/2.0.6
##module load tophat/2.0.13
##module load bedtools/v2.25.0

rsemDir=$(dirname `which rsem-prepare-reference`)

echo begin at: `date` `hostname`

### run hisat2
$sDir/../rna/aln.hisat2.sh -t $optT $outDir/$sampleID $fq1 $fq2 $sampleID "L001" "001"
samtools view -b -f 0x4 $outDir/$sampleID/$sampleID.hisat.hit.sort.bam -o $outDir/$sampleID/$sampleID.hisat.unmapped.bam
### trapes
### some problem with samtools/1.2.0
module load samtools/1.3.1
python /Share/BP/zhenglt/01.bin/RepSeq/TRAPeS/trapes.py \
    -rsem $rsemDir \
    -genome hg19 \
    -lowQ \
    -path $outDir \
    -sumF $outDir/$sampleID.trapes.TCR.sum \
    -output $sampleID.trapes.TCR.out \
    -bam $sampleID.hisat.hit.sort.bam \
    -unmapped $sampleID.hisat.unmapped.bam
###    -bowtie2 \
###    -samtools \


echo end at: `date`
