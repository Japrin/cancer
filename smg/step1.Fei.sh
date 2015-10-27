#!/bin/bash
source /PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh

echo begin at: `date` $$

wDir=/BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/smg
cd $wDir
outDir=coverage.Fei
mkdir -p $outDir/gene_coverage
echo REF: $REF

genome music bmr calc-covg \
	--roi-file roi.00.txt \
	--reference-sequence  $REF \
	--bam-list $wDir/bam.S12.Fei.list \
	--output-dir $outDir \
	--gene-covg-dir $outDir/gene_coverage \
	--normal-min-depth 8 \
	--tumor-min-depth 8 \
	--min-mapq 1 
	#--cmd-list-file step1.job
	#--cmd-prefix "qsub"
	#--cmd-prefix "qsub -l vf=1g -cwd -V"

echo end at: `date`
