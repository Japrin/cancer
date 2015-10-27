#!/bin/bash -eu

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
method="general"
optA=""
optB="INDEL"
gender="M"

while getopts c:m:a:b:f:g: opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid configure file ($OPTARG), default will be used"
		fi
		;;
	m)
		method="$OPTARG"
		;;
	a)
		optA="$OPTARG"
		;;
	b)
		optB="$OPTARG"
		;;
	g)
		if [ -f $OPTARG ]
		then
			source "$OPTARG"
		else
			gender=$OPTARG
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-m method, default 'general'] [-a 'somatic' or ''] [-g gender|gender file] <infile> <ID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile] [-m method, default 'general'] [-a 'somatic' or ''] [-g gender|gender file] <infile> <ID>"
	exit 1
fi

echo begin at: `date`

source $iniFile 

infile=$1
sampleID=$2
outDir=`dirname $infile`

echo ">>> annotate ..."

mkdir -p $outDir

cancerGeneDir="/Share/BP/zhenglt/00.database/cancer_gene"

#### read from stdin, and ouput to stdout
function addGeneAnn
{
	perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName OMIM_DAVID \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName CancerGene \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName BertVogelstein125 -annFile $cancerGeneDir/compile/Gene.manual.BertVogelstein125 \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName Predisposition -annFile $cancerGeneDir/cancer_gene_predisposition.slim.txt \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName DriverCNA -annFile $cancerGeneDir/cancer_gene_cnv.slim.txt \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName Rearrangement -annFile $cancerGeneDir/cancer_gene_rearrangment.slim.txt \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName GO_BP \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName GO_CC \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName GO_MF \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName KEGG_PATHWAY \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName PANTHER_PATHWAY \
	| perl $PIPELINE/cancer/huyy/addAnn_hyy.pl -annName REACTOME_PATHWAY
}


if [[ "$infile" =~ ".vcf" || "$infile" =~ ".vcf.gz" ]]; then

	## annotate with ANNOVAR
	annovarFile=${infile/.vcf*/.annovar}
	addAnnFile=${infile/.vcf*/.annovar.hg19_multianno.txt}
	vcfFile=${infile/.vcf*/.reformated.vcf}
	## remove Y chromosomes in FEMALE cases And normolize
	if [[ "$gender" =~ "F" ]]; then
		bcftools filter -t ^Y,chrY $infile -i 'FILTER=="PASS"||FILTER=="."' \
			| bcftools norm -m-both \
			| bcftools norm -f $REF \
			| convert2annovar.pl --withzyg --includeinfo -format vcf4old - \
			> $annovarFile
	else
		bcftools filter $infile -i 'FILTER=="PASS"||FILTER=="."' \
			| bcftools norm -m-both \
			| bcftools norm -f $REF \
			| convert2annovar.pl --withzyg --includeinfo -format vcf4old - \
			> $annovarFile
	fi

	if [ ! -f "$addAnnFile" ];then
		table_annovar.pl $annovarFile $HumanDB -buildver hg19 -otherinfo -remove -nastring . \
			-protocol knownGene,cytoband,genomicSuperDups,phastConsElements46way,wgRna,targetScanS,tfbsConsSites,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,avsnp142,ljb26_all,gerp++elem,gerp++gt2,cosmic70 \
			-operation g,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f
	fi

	## refomrat the ouput table_variants file, mainly for the otherinfo field
	bcftools view -h $infile | awk '/^##/' > $infile.header.tmp
	perl  $PIPELINE/cancer/var/reformat_annovar.pl -id $sampleID $addAnnFile \
		| perl $PIPELINE/cancer/var/tab2vcf.pl -id $sampleID -f $infile.header.tmp \
		| bcftools norm -m+both \
		> $vcfFile
	rm $infile.header.tmp	
	## re-annoate indel (not nessesery if left-normed
	######if [[ "$optB" =~ "INDEL" ]] || [[ "$optB" =~ "indel" ]]; then
	######	perl $PIPELINE/var/samtools-dbsnp.pl -m INDEL -i $vcfFile -o $vcfFile.tmp1
	######	perl $PIPELINE/var/samtools-1000G.pl  -i $vcfFile.tmp1 -o $vcfFile.tmp2
	######	perl $PIPELINE/var/samtools-cosmic.pl -i $vcfFile.tmp2 -o $vcfFile.tmp3
	######	mv $vcfFile.tmp3 $vcfFile
	######	rm $vcfFile.tmp1 $vcfFile.tmp2
	######fi

	## somatic filter
	if [[ "$optA" =~ "somatic" ]]; then
		echo "Somatic filter"
		perl -i -F"\t" -ane 'if(/^#/ || ($F[7]!~/snp142/ && $F[2] eq".") ||  ($F[7]=~/cosmic/ && ($F[7]!~/1000g2014oct_all=(.+?);/ || ($F[7]=~/1000g2014oct_all=(.+?);/ && $1<0.01)) ) ){print}' $vcfFile
	else
		echo "No somatic filter"
	fi
	
	vcfgz=${vcfFile/.vcf/.vcf.gz}
	bgzip -f $vcfFile
	tabix -f -p vcf $vcfgz
	rm -f $annovarFile

	## stat
	vcf.func.stat.pl -f $vcfgz > $vcfgz.stat.txt

elif [[ "$infile" =~ ".bed" ]]; then
	echo annotate bed file: $infile
	awk -F"\t" -v OFS="\t" '{print $1,$2-1,$3,0,0,"Name="$4}' $infile > $infile.avinput
	table_annovar.pl $infile.avinput $HumanDB -buildver hg19 -otherinfo -remove -nastring . \
		-protocol knownGene,evofold,wgRna,targetScanS,phastConsElements46way,genomicSuperDups,tfbsConsSites,cytoband,gff3 \
		-operation g,r,r,r,r,r,r,r,r \
		--gff3dbfile hg19_rmsk.gff
	echo annova done, see result: $infile.avinput.hg19_multianno.txt
	rm $infile.avinput

elif [[ "$infile" =~ ".gff" ]]; then

	###"SVID" "SVType must be in column9;
	awk -F"\t" -v OFS="\t" '{print $1,$4,$5,"0","0",$9;}' $infile > $infile.avinput

	table_annovar.pl $infile.avinput $HumanDB -buildver hg19 -otherinfo -remove -nastring . \
		-protocol knownGene,evofold,wgRna,targetScanS,phastConsElements46way,genomicSuperDups,tfbsConsSites,cytoband,gff3 \
		-operation g,r,r,r,r,r,r,r,r \
		--gff3dbfile hg19_rmsk.gff

	perl $PIPELINE/cancer/var/multianno.reformat.pl $infile.avinput.hg19_multianno.txt | sed '1,1s/gff3/Repeat/' | sed '1,1s/.knownGene//g' > $infile.ann.xls
	perl $PIPELINE/cancer/var/var_sv_putative.fusion.gene.pl $infile.ann.xls > $infile.ann.fusionGene.xls

	## all gene 
	perl $PIPELINE/cancer/var/cnv.geneInfo.pl $infile.ann.xls \
		| addGeneAnn \
		> $infile.ann.geneInfo.xls

	## only fusion genes
	cat $infile.ann.fusionGene.xls \
		| addGeneAnn \
		> $infile.ann.fusionGene.full.xls

        #$infile.avinput.hg19_multianno.txt 
	rm -f $infile.avinput $infile.avinput.*.invalid_input 2>/dev/null
fi

echo "*** Finished annotating variants ***"
echo end at: `date`
