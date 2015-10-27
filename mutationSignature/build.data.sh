

awk -v OFS="\t" '{print $3,$5,$6,$13,$2,$4}' /WPS/GR/zhengliangtao/01bin/annovar/annovar_2013Feb25/humandb_hg19/hg19_refGene.txt | sort -k 1,1 -k 2g,2 -k 3g,3 > hg19_refGene.transcript.txt
