#!/bin/bash

for i in {1..22} X
do
	echo $i
	outDir=/WPS/GR/zhengliangtao/01bin/novoHumanReseq_v1.0.1/app/novoHumanReseq/bin/cancer/mutationSignature/context/chr$i
	mkdir -p $outDir
	python /WPS/GR/niuzexiong/program/fa_deal/sum_N_combinat_depth.py \
		-d 10 \
		-r /WPS/GR/zhengliangtao/work/BT/BT_21_ganai/CapData/TruSeq-Custom-Enrichment-Trial-Kit-Regions-bed-file.slim.bed \
		-o $outDir \
		-a /WPS/GR/zhengliangtao/00database/human/GATK_bundle/1.5/hg19/byChr/chr$i.fa \
		-i /WPS/GR/zhengliangtao/work/zigongneimo_NH130032/analysis_20130220/aln/TZigongneimo/all/byChr/chr$i/TZigongneimo.chr$i.recal.bam \
		-q 13 -Q 1

done
