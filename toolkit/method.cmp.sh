#!/bin/bash

outDir=somatic.snv.cmp
while read id normalBam tumorBam
do
	echo $id...
	mkdir -p $outDir/$id
	venn15.pl $outDir/$id -i=0,1 -n1 muTect -n2 strelka -n3 varScan -n4 samtools \
		../muTect/out/$id/$id.muTect.call.filter.reformated.vcf.gz \
		../strelka/out/$id/results/$id.strelka.somatic.snvs.filter.reformated.vcf.gz \
		../varScan/out/$id/$id.varScan.call.Somatic.reformated.vcf.gz \
		../samtools/out/$id/$id.samtools.somatic.snv.fpfilter.reformated.vcf.gz \
		-title $id.somatic.snv
	java -Djava.awt.headless=true -jar /PROJ/HEALTH/zhengliangtao/03toolkit/svg_kit/batik-1.7/batik-rasterizer.jar \
		-m image/png \
		$outDir/$id/$id.somatic.snv.svg
	cp $outDir/$id/$id.somatic.snv.png $outDir
done<../bam.list


outDir=somatic.indel.cmp
while read id normalBam tumorBam
do
	echo $id...
	mkdir -p $outDir/$id
	venn7.pl $outDir/$id -i=0,1 -n1 strelka -n2 varScan -n3 samtools \
		../strelka/out/$id/results/$id.strelka.somatic.indels.filter.reformated.vcf.gz \
		../varScan/out/$id/$id.varScan.indel.Somatic.hc.reformated.vcf.gz \
		../samtools/out/$id/$id.samtools.somatic.indel.reformated.vcf.gz \
		-title $id.somatic.indel
	java -Djava.awt.headless=true -jar /PROJ/HEALTH/zhengliangtao/03toolkit/svg_kit/batik-1.7/batik-rasterizer.jar \
		-m image/png \
		$outDir/$id/$id.somatic.indel.svg
	cp $outDir/$id/$id.somatic.indel.png $outDir
done<../bam.list
