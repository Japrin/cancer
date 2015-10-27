#!/bin/bash


source /PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh

cd /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/smg
echo begin at: `date`


outDir=pathway.non-normal.G1
mkdir -p $outDir

while read ID f
do
	genome music path-scan --pathway-file $f --gene-covg-dir coverage.non-normal.G1/gene_covgs --bam-list bam.S23.non-normal.list --maf-file /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/strelka/somatic_mutation.non-normal.S23.maf  --output-file $outDir/pathway.$ID --min-mut-genes-per-path 2 --bmr 4.0e-6

	genome music path-scan --pathway-file $f --gene-covg-dir coverage.non-normal.G1/gene_covgs --bam-list bam.S23.non-normal.list --maf-file /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/strelka/somatic_mutation.non-normal.S23.maf  --output-file $outDir/pathway.$ID.exclude.TP53 --min-mut-genes-per-path 2 --bmr 4.0e-6 --genes-to-ignore TP53


done<path.list

echo end at: `date`
