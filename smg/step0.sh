cat \
	/PUBLIC/database/HEALTH/annovar/humandb_b37/ucscGene/hg19_ucscGene.CDS.b37.bed \
	/PUBLIC/database/HEALTH/annovar/humandb_b37/ucscGene/hg19_ucscGene.Splicing.b37.bed \
	| awk -F"\t" -v OFS="\t" '$3-$2>0 && $1!~/[^0-9XY]/ {print $1,$2,$3,$4}' \
	| sort -k 1,1 -k 2g,2 -k 3g,3 \
	| mergeBed -i stdin -nms \
	| perl -F"\t" -ane 'chomp @F;@aa=split /,/,$F[3];$F[1]++;%bb=();map {$bb{$_}++} @aa;print join("\t",@F[0,1,2])."\t".join(",",(keys %bb)[0])."\n";' \
	> roi.00.txt
