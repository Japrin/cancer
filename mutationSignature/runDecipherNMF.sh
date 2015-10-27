
##
#bgzip -cd /WPS/BP/zhenglt/work/novogene/breastCancer/release/somatic_snv/2711B-vs-109PT-2.genome.final.vcf.gz | /WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/mutationSignature/snv_context_96.pl -f /WPS/BP/zhenglt/00.database/broad/bundle/2.8/hg19/hg19.fa -s 109PT-2 > 109PT-2.context96.txt
#
#bgzip -cd /WPS/BP/zhenglt/work/novogene/breastCancer/release/somatic_snv/2711B-vs-189NT.genome.final.vcf.gz   | /WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/mutationSignature/snv_context_96.pl -f /WPS/BP/zhenglt/00.database/broad/bundle/2.8/hg19/hg19.fa -s 189NT   > 189NT.context96.txt
#
#bgzip -cd /WPS/BP/zhenglt/work/novogene/breastCancer/release/somatic_snv/2711B-vs-332NT-1.genome.final.vcf.gz | /WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/mutationSignature/snv_context_96.pl -f /WPS/BP/zhenglt/00.database/broad/bundle/2.8/hg19/hg19.fa -s 332NT-1 > 332NT-1.context96.txt


##  n

#ls data_Nik-ZainalEtAl/bySample/*.context96 \
#	| xargs paste 109PT-2.context96.txt 332NT-1.context96.txt \
#	| awk -F"\t" -v OFS="\t" 'BEGIN{ printf "#1.2\n96\t23\n";} { printf $1"\t"$1; for(i=2;i<=NF;i+=2){printf "\t"$i} printf "\n"}' \
#	| sed 's/context/ID1/' | sed 's/context/ID2/' \
#	> S23.context96.gct

#sed '1,2d' S23.context96.gct \
#	| perl -F"\t" -ane 'chomp @F;if($F[0]=~/.\((.+?)\)./){$F[0]=$1} if($F[1]=~/(.)\((.)>.\)(.)/){$F[1]="$1$2$3"}  print join("\t",@F)."\n";' \
#	| sed 's/ID1/types/' | sed 's/ID2/subtypes/' \
#	> S23.context96.txt
odir=matDir.nfact5
mkdir -p  $odir
cut -f1  S23.context96.txt > $odir/types.txt
cut -f2  S23.context96.txt > $odir/subtypes.txt
head -1 S23.context96.txt | transpose.sh | sed '1,2d' | awk 'BEGIN{print "sampleName"}{print}' > $odir/sampleNames.txt
sed '1,1d' S23.context96.txt | cut -f 1-2 --complement > $odir/originalGenomes.txt

mkdir -p temp
mkdir -p out
/WPS/BP/zhenglt/01.bin/matlab/R2013a/mybuild/bin/matlab -nosplash -nodisplay -r "addpath('/WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/mutationSignature'); runDecipherNMF('$odir',5);quit"
(
	sed '1,1d' $odir/sampleNames.txt | transpose.sh | xargs echo "Name" "Description" | tr ' ' '\t'
	awk -F"\t" -v OFS="\t" '{print "Sig"NR,"Sig"NR,$0 }' $odir/S.H
) > $odir/S.H.gct

(
	printf "types\tsubtypes\tSig1\tSig2\tSig3\tSig4\tSig5\n"
	paste <(cut -f 1,2 S23.context96.gct | sed '1,3d') $odir/S.W
) > $odir/S.W.gct

/WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/mutationSignature/NMF.plot.R $odir/S.W.gct $odir/S.W.png $odir/S.H.gct $odir/S.H.png
