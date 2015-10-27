#!/bin/bash

source /PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh

cd /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/smg

echo begin at: `date`
genome music bmr calc-bmr \
	--bam-list /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/smg/bam.S23.non-normal.list \
	--reference-sequence $REF \
	--roi-file /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/smg/roi.00.txt \
	--output-dir /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/smg/coverage.non-normal.G1 \
	--bmr-groups 1 \
	--maf-file /BJPROJ/GR/HUMAN/sci/WES.sumin.exome/mycall/strelka/somatic_mutation.non-normal.S23.maf \
	--bmr-output mr.non-normal.G1/bmr.txt \
	--gene-mr-file mr.non-normal.G1/gene-mr.txt 

echo end at: `date`
