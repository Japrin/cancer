#!/bin/bash
echo "*** report_personel.sh ***"

iniFile="/WPS/GR/zhengliangtao/02pipeline/novo.med.seq/cancer/parameter/init_human.sh"

while getopts c: opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid reference file ($OPTARG), default will be used"
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] <sampleID> <inBam> <snpvcf.gz> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))



if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] <sampleID> <inBam> <snpvcf.gz> <outDir>"
	exit 1
fi

sampleID=$1
in_bam=$2
snpVCF=$3
outDir=$4

echo begin at: `date` 
source $iniFile
###----------------------- GWAS ----------------------------###

targetRegion="/WPS/GR/zhengliangtao/01bin/annovar/annovar_2013May09/humandb_b37/b37_gwasCatalog.site.bed"

mkdir -p $outDir

varDir="$outDir/gwas"
mkdir -p $varDir
raw_snp=$varDir/$sampleID.GATK.gwas.genotype.raw.vcf
flt_snp=$varDir/$sampleID.GATK.gwas.genotype.flt.vcf

samtools index $in_bam
## snp calling
java -Djava.io.tmpdir=$varDir -Xmx5G -jar $gatkJAR -T UnifiedGenotyper \
	-out_mode EMIT_ALL_CONFIDENT_SITES \
	-L $targetRegion \
	-R $refData \
	-I $in_bam \
	-o $raw_snp \
	-stand_call_conf 30 \
	-stand_emit_conf 10.0 \
	-A QualByDepth \
	-A FisherStrand \
	-A AlleleBalance \
	-A DepthOfCoverage \
	-A MappingQualityZero \
	-A TandemRepeatAnnotator \
	-baq CALCULATE_AS_NECESSARY \
	-glm SNP \
	-et NO_ET \
	-K $gatkKey
echo finish snp calling at: `date`

## basic snp filtering
filter="QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
java -Djava.io.tmpdir=$varDir -Xmx5G -jar $gatkJAR -T VariantFiltration \
	-L $targetRegion \
	-R $refData \
	-o $flt_snp \
	-V $raw_snp \
	--filterExpression "$filter" \
	--filterName "StandardFilter" \
	-et NO_ET \
	-K $gatkKey
echo "finish variantFiltration(snp)" at: `date`

## gwas search

inVCF="$flt_snp"
perl -F"\t" -ane '$line=$_;chomp @F;if(/^#/){next};printf "$F[0]\t$F[1]\t$F[1]\t$F[3]\t%s\t%s",$F[4] eq "."?$F[3]:$F[4],$line;' $inVCF > $inVCF.annovar
annotate_variation.pl -buildver hg19 -regionanno -dbtype gwasCatalog $inVCF.annovar $HumanDB
awk -F "\t" '$3!~/_/' $inVCF.annovar.hg19_gwasCatalog \
	| perl -F"\t" -ane '/DP=(.+?);/;if($1>=4 && $F[12]>=10){printf "$_";}' \
	| perl -F"\t" -ane 'chomp @F;$F[1]=~s/:(\d{4,}),/:\1\|/g;printf "%s\n",join("\t",@F);' \
	> $inVCF.annovar.hg19_gwasCatalog.confident.txt
perl -F"\t" -ane 'if(/:(rs\d+)-(\w):/ && $2 eq $F[6]){print "$_"};' $inVCF.annovar.hg19_gwasCatalog.confident.txt > $inVCF.annovar.hg19_gwasCatalog.confident.risk.txt
perl -ane '/Name=(.+?):/;print "$1\t$_";' $inVCF.annovar.hg19_gwasCatalog.confident.risk.txt \
	| sort -k 1,1 | cut -f 1 | uniq -c | sort -k 1gr,1 \
	> $inVCF.annovar.hg19_gwasCatalog.confident.risk.trait.list
awk '/diabet|cancer/' $inVCF.annovar.hg19_gwasCatalog.confident.risk.trait.list \
	| perl -ane 's/^\s+\d+\s//;print "$_"' \
	> $inVCF.annovar.hg19_gwasCatalog.confident.risk.trait.list.example
perl /WPS/GR/zhengliangtao/01bin/novoHumanReseq_v1.0.1/app/novoHumanReseq/bin/gwas.example.pl \
	-trait $inVCF.annovar.hg19_gwasCatalog.confident.risk.trait.list.example \
	$inVCF.annovar.hg19_gwasCatalog.confident.risk.txt \
	> $inVCF.annovar.hg19_gwasCatalog.confident.risk.example.xls

ofile="$varDir/$sampleID.gwas.summary.txt"
printf "### summary of gwas searching ###\n" > $ofile
t=`cat $inVCF.annovar.hg19_gwasCatalog.confident.txt | wc -l`
printf "high confidence genotype sites: $t\n" >> $ofile 
t=`cat $inVCF.annovar.hg19_gwasCatalog.confident.risk.txt | wc -l`
printf "high confidence sites with risk allele: $t\n" >> $ofile 
t=`cat $inVCF.annovar.hg19_gwasCatalog.confident.risk.trait.list | wc -l`
printf "associated phonotypes: $t\n" >> $ofile 
t=`sed '1,1d' $inVCF.annovar.hg19_gwasCatalog.confident.risk.example.xls | cut -f 2,3 | sort | uniq | wc -l`
printf "associated with diabet or cancer: $t\n" >> $ofile 

echo finish gwas searching at `date`

###----------------------- COSMIC ----------------------------###
tBed="/WPS/GR/zhengliangtao/00database/cancer/COSMIC/data_export/CosmicCompleteExport_v60_190712.slim.SNV.bed"
mkdir -p $outDir/cosmic
tOut=$outDir/cosmic/$sampleID.GATK.snp.COSMIC.out
tXls=$outDir/cosmic/$sampleID.GATK.snp.COSMIC.xls

bgzip -cd $snpVCF | awk '!/^#/ && $7=="PASS" ' \
	| intersectBed -a stdin -b $tBed -wa -wb \
	> $tOut

perl -ane  'if(!/Func=(exonic|splicing)/){ next; } if(/ExonicFunc=synonymous_SNV/){ next; } if(/ESP.+?_ALL=(.+?);/) { $esp=$1; } if(/1000g.+?_ALL=(.+?);/) { $kg=$1; } if($esp<0.01 || $kg<0.01){ print }' $tOut \
	> $tOut.example.txt
#chr1	26515956	rs1045105	C	A	0.118423	0.10	CNKSR1	NM_006314	c.C2059A	p.H687N	+	Confirmed somatic variant	CNKSR1	PD2125a	c.2059C>T	p.H687Y	Substitution - Missense	het	20054297
printf "chr\tpos\tdbSNP ID\tref\talt\tESP\t1KG\tGene\tTranscript\tcDNA\taa change\tstrand\tstatus\tGene (cosmic)\tsample (cosmic)\tcDNA (cosmic)\taa change (cosmic)\ttype\tzygosity (cosmic)\tpubmed ID\n" > $tXls
perl -F"\t" -ane 'chomp @F;$nF=$#F;($esp)=/ESP.+?_ALL=(.+?);/;($KG)=/1000g.+?_ALL=(.+?);/;if(/Gene=(.+?);/){$gene=$1}else{$gene="na"};;if(/AAChange=(.+?);/){$AAChange=$1;$AAChange=~s/:/\t/g;printf "%s\t%s\t%s\t$gene\t$AAChange\t%s\n",join("\t",@F[0..4]),$esp?$esp:"na",$KG?$KG:"na",join("\t",@F[13..$nF]);}else{printf "%s\t%s\t%s\t$gene\tna\tna\tna\t%s\n",join("\t",@F[0..4]),$esp?$esp:"na",$KG?$KG:"na",join("\t",@F[13..$nF]);}' $tOut >> $tXls
perl -F"\t" -ane 'chomp @F;$nF=$#F;($esp)=/ESP.+?_ALL=(.+?);/;($KG)=/1000g.+?_ALL=(.+?);/;if(/Gene=(.+?);/){$gene=$1}else{$gene="na"};;if(/AAChange=(.+?);/){$AAChange=$1;$AAChange=~s/:/\t/g;printf "%s\t%s\t%s\t$gene\t$AAChange\t%s\n",join("\t",@F[0..4]),$esp?$esp:"na",$KG?$KG:"na",join("\t",@F[13..$nF]);}else{printf "%s\t%s\t%s\t$gene\tna\tna\tna\t%s\n",join("\t",@F[0..4]),$esp?$esp:"na",$KG?$KG:"na",join("\t",@F[13..$nF]);}' $tOut.example.txt >> $tOut.example.xls

ofile="$outDir/cosmic/$sampleID.cosmic.summary.txt"
printf "### summary of cosmic searching ###\n" > $ofile
t=`cut -f 1,2 $tOut | sort | uniq | wc -l`
printf "snv recorded in COSMIC: $t\n" >> $ofile
t=`cut -f 1,2 $tOut.example.txt | sort | uniq | wc -l`
printf "low frequency snv recorded in COSMIC: $t\n" >> $ofile

###----------------------- PharmGKB ----------------------------###
mkdir -p $outDir/pharmGKB
bgzip -cd $snpVCF \
	| perl -ane 'BEGIN{open IN,"/WPS/GR/zhengliangtao/00database/human/PharmGKB/clinicalAnnotations.txt";while(<IN>){chomp;@F=split /\t/;$d{$F[0]}=$_;} } foreach $k (keys %d) {if(/\t$k\t/){printf "%s\t%s",$d{$k},$_;delete $d{$k}; }} ' \
	> $outDir/pharmGKB/$sampleID.pharmGKB.txt

ofile="$outDir/pharmGKB/$sampleID.pharmGKB.summary.txt"
printf "### summary of pharmGKB searching ###\n" > $ofile
t=`cat $outDir/pharmGKB/$sampleID.pharmGKB.txt | wc -l`
printf "snp with PharmGKB record: $t\n" >> $ofile

###----------------------- END ----------------------------###

echo end at `date`
echo "*** finish report_personel.sh ***"


