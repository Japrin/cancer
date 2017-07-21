#!/bin/bash

if [ $# -lt 4 ]
then 
	echo "Usage: $0 <outDir> <sampleID> <normalBam> <tumorBam>"
	exit 1
fi

outDir=$1
sampleID=$2
normalBam=$3
tumorBam=$4

optT=4

mkdir -p $outDir

. /usr/share/Modules/init/bash
export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
module load samtools/0.1.19
export PSHOME=/Share/BP/zhenglt/01.bin/hla/POLYSOLVER/polysolver

echo begin at: `date`
mkdir -p $outDir
source $PSHOME/scripts/config.bash
####$PSHOME/scripts/shell_call_hla_type </path/to/bam> <race> <includeFreq> <build> <format> <insertCalc> </path/to/output_directory>
$PSHOME/scripts/shell_call_hla_type $normalBam Unknown 1 hg19 STDFQ 0 $outDir >$outDir/$sampleID.polysolver.log 2>&1
echo "finish hlatyping at: " `date`

####$PSHOME/scripts/shell_call_hla_mutations_from_type </path/to/normal_bam> </path/to/tumor_bam> </path/to/winners.hla.txt> <build> <format> </path/to/output_directory>
$PSHOME/scripts/shell_call_hla_mutations_from_type $normalBam $tumorBam $outDir/winners.hla.txt hg19 STDFQ $outDir >>$outDir/$sampleID.polysolver.log 2>&1
echo "finish hla mutation calling at: " `date`

####$PSHOME/scripts/shell_annotate_hla_mutations <prefix_to_use> </path/to/directory_with_mutation_detection_output>
$PSHOME/scripts/shell_annotate_hla_mutations indiv $outDir >>$outDir/$sampleID.polysolver.log 2>&1
echo "finish hla mutation annotating at: " `date`

echo end at: `date`
