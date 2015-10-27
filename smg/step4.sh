#!/bin/bash


source /PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh

cd /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/smg
echo begin at: `date`


outDir=pathway
mkdir -p $outDir

while read ID f
do
	genome music path-scan --pathway-file $f --gene-covg-dir coverage/gene_covgs --bam-list bam.S12.Ca.list --maf-file /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/strelka/somatic_mutation.Ca.S12.maf  --output-file $outDir/pathway.$ID --min-mut-genes-per-path 2 --bmr 4.0e-6

	genome music path-scan --pathway-file $f --gene-covg-dir coverage/gene_covgs --bam-list bam.S12.Ca.list --maf-file /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/strelka/somatic_mutation.Ca.S12.maf  --output-file $outDir/pathway.$ID.exclude.TP53 --min-mut-genes-per-path 2 --bmr 4.0e-6 --genes-to-ignore TP53


done<path.list

echo end at: `date`
