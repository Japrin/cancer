

## symbol, mRNA accesion, protein accesion, entrez id
perl -F"\t" -ane 'chomp @F;if($F[2]=~/^NR_/){next} @a=split /; /,$F[0]; foreach (@a) { printf "%s\n",join("\t",$_,@F[2,3,6]); }' \
	/WPS/GR/zhengliangtao/01bin/annovar/annovar_2013Feb25/humandb_hg19/hg19_refLink.txt \
	| sort -t $'\t' -k 3,3 \
	> ./hg19.idmapping.1

gzip -cd /WPS/GR/zhengliangtao/00database/uniprot/knowledgebase/idmapping/HUMAN_9606_idmapping_selected.tab.gz \
	| perl -F"\t" -ane 'chomp @F;if($F[3] eq ""){next} @a=split /; /,$F[3]; foreach (@a){ s/\.\d+$//; printf "%s\n",join("\t",$F[0],$_) }' \
	| sort -t $'\t' -k 2,2 \
	> ./NP.UniProt.mapping

## NP_001108225    ENG     NM_001114753    2022    P17813
## protein accesion, symbol, mRNA accesion, Entrez ID, UniProtKB-AC
join -t $'\t' -1 3 -2 2 hg19.idmapping.1 NP.UniProt.mapping | sort -k 3,3 -k 5,5  > hg19.idmapping.2
