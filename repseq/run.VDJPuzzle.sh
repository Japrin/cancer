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

mkdir -p $outDir/input/$sampleID
#sampleID=SAMPLEID
#inDir=/Share/BP/zhenglt/01.bin/RepSeq/VDJPuzzle/Example
#outDir=/Share/BP/zhenglt/01.bin/RepSeq/VDJPuzzle/test

###source $iniFile

export MODULESHOME=/usr/share/Modules
. /usr/share/Modules/init/bash
export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
module load novomedSeq/1.0
module load igblast/1.4.0
module load java/1.7.0_79
module load trinity/2.0.6
module load tophat/2.0.13
module load bowtie2/2.2.3
module load samtools/0.1.19
module load bedtools/v2.25.0

VDJPuzzleBin=/Share/BP/zhenglt/01.bin/RepSeq/VDJPuzzle/VDJPuzzle.sh

echo begin at: `date` `hostname`


ln -s $fq1 $outDir/input/$sampleID/${sampleID}_L001_R1_001.fastq.gz
ln -s $fq2 $outDir/input/$sampleID/${sampleID}_L001_R2_001.fastq.gz
cd $outDir
$VDJPuzzleBin $outDir/input > $outDir/$sampleID.VDJPuzzle.log 2>&1

echo end at: `date`
