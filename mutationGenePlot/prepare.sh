
infile="/WPS/GR/niuzexiong/project/cancer-liver-exome/02.analysis/130406hiseq/whole_analy0509/316/somatic/varScan/316.varScan.somatic.snv.exome.final.txt"
outfile="S316.snv.final.for.MutationGenePlot.txt"
tmpfile="./tmp"
mappingfile="/WPS/GR/zhengliangtao/01bin/novoHumanReseq_v1.0.1/app/novoHumanReseq/bin/cancer/mutationGenePlot/data/hg19.idmapping.2"

sed '1,1d' $infile \
	|  perl -F"\t" -ane 'chomp @F; if($F[10]=~/NA|UNKNOWN/){next}; $F[8]=~s/[\[\]]//g; @a=split /, /,$F[8]; foreach (@a){ s/:/\t/g; printf "$F[0]\t$F[1]\t$_\n"; }  ' \
	| sort -t $'\t' -k 4,4 \
	> $tmpfile

join -t $'\t' -1 4 -2 3 $tmpfile $mappingfile \
	| sort -t $'\t' -k 4,4 -k 11,11 -k 2,2 -k 3g,3 \
	| awk -F"\t" -v OFS="\t" '{print $0,$2,$3,$11}' |  uniq -f 11 | cut -f 1-11 \
	> $outfile
