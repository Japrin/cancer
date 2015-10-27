


#./snv_spectrum.pl -p \
#	/WPS/GR/zhengliangtao/work/zigongneimo_NH130032/analysis_20130220/somatic/Zigongneimo/varScan/all.varScan.somatic.snv.genome.final.vcf
#
#./addStrand.pl \
#	-l hg19_refGene.transcript.txt \
#	/WPS/GR/zhengliangtao/work/zigongneimo_NH130032/analysis_20130220/somatic/Zigongneimo/varScan/all.varScan.somatic.snv.genome.final.vcf \
#	| ./strandbias_spectrum.pl -p


./surroundingSequence.pl /WPS/GR/zhengliangtao/00database/human/GATK_bundle/1.5/hg19/hg19.fa \
	/WPS/GR/zhengliangtao/work/zigongneimo_NH130032/analysis_20130220/somatic/Zigongneimo/varScan/all.varScan.somatic.snv.genome.final.vcf \
	-s 5 -o test.out -p -triplet
