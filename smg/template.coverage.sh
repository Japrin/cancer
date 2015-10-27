#!/bin/bash

sampleID=$1
normalBam=$2
tumorBam=$3

source /PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh
wDir=/BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/smg
cd $wDir
outDir=coverage
mkdir -p $outDir/gene_coverage
echo REF: $REF

echo begin at: `date`

mkdir -p $outDir/roi_covgs

genome music bmr calc-covg-helper \
	--normal-tumor-bam-pair "$sampleID	$normalBam	$tumorBam" \
	--roi-file "$wDir/roi.00.txt" \
	--reference-sequence "$REF" \
	--output-file "$outDir/roi_covgs/$sampleID.covg" \
	--normal-min-depth 8 \
	--tumor-min-depth 8 \
	--min-mapq 1

echo end at: `date`
