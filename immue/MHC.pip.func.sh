#!/bin/bash

function vcf2annovar
{
    infile=$1
    outDir=$2
    outfile=$outDir/`basename $infile`
    if [[ "$infile" =~ ".vcf.gz" ]]; then
        #echo vcf.gz file...
        outfile=${outfile/.vcf.gz/.annovar}
        bgzip -cd $infile | convert2annovar.pl --withzyg --includeinfo -format vcf4old - > $outfile
    elif [[ "$infile" =~ ".vcf" ]]; then
        #echo vcf file...
        outfile=${outfile/.vcf/.annovar}
        convert2annovar.pl --withzyg --includeinfo -format vcf4old $infile > $outfile
    fi
    echo $outfile
}

function run
{
	echo begin at: `date`
	. /usr/share/Modules/init/bash
	export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
	module load novomedSeq/1.0

	humanDB=/DBS/DB_temp/zhangLab/annovar/humandb_b37
	geneModel=knownGene
	mRNAGenePred=$humanDB/hg19_knownGene.txt
	mRNAFa=$humanDB/hg19_knownGeneMrna.fa

	ifile=$1
	outDir=$2
	sampleID=$3
	chr=$4
	kmer=$5

	prefix=$outDir/`basename $ifile`.$chr.k$kmer

	mkdir -p $outDir
	if [ "$chr" == "all" ];then
		ln -s $ifile $prefix
	else
		awk '$1=="'$chr'"' $ifile > $prefix
	fi
	annotate_variation.pl --geneanno --buildver hg19 --dbtype $geneModel --separate $prefix $humanDB

	#/WPS1/zhenglt/work/MHC/TCGA/bin/annovar.variant.codingSequence.pl --verbose \
	/WPS1/zhenglt/work/MHC/TCGA/bin/annovar.variant.codingSequence.pl \
		-i \
		--kmer $kmer \
		--peptide $prefix.exonic_variant_function.codingSeq.peptide.fa \
		--protein $prefix.exonic_variant_function.codingSeq.protein.fa \
		--sample $sampleID \
		$prefix.exonic_variant_function \
		$mRNAGenePred \
		$mRNAFa \
		> $prefix.exonic_variant_function.codingSeq.out \
		2> $prefix.exonic_variant_function.codingSeq.err

	echo end at: `date`
}

function runIEDB
{
	id=$1
	pfile=$2
	outDir=$3
	hla=$4
	plen=$5
	
	echo .......... begin at: `date` "($id $hla $plen)" ..........
	. /usr/share/Modules/init/bash
	export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
	module load novomedSeq/1.0
	export NMHOME=/Share/BP/zhenglt/01.bin/netMHC/netMHC-3.4
	export TMPDIR="$outDir/$id.$hla.$plen/tmp"
	mkdir -p $TMPDIR
	chmod 1777 $TMPDIR
	python /Share/BP/zhenglt/01.bin/IEDB/mhc_i/src/predict_binding.py IEDB_recommended "$hla" $plen "$pfile" > "$outDir/IEDB.$id.$hla.$plen.out"
	rm -r $outDir/$id.$hla.$plen
	echo .......... end at: `date`  "($id $hla $plen)" ..........

}
function runNetMHC
{
	id=$1
	pfile=$2
	outDir=$3
	hla_allele="HLA-A02:01"
	plen=9
	
	echo begin at: `date` hla: $hla_allele len: $plen
	. /usr/share/Modules/init/bash
	export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
	module load novomedSeq/1.0

	##hla_list_file="/WPS1/zhenglt/work/MHC/TCGA/netMHC.HLA.list"
	##/Share/BP/zhenglt/01.bin/netMHC/netMHC-3.4/netMHC -A | awk '/^HLA/' > $hla_list_file

	export NMHOME=/Share/BP/zhenglt/01.bin/netMHC/netMHC-3.4
	export TMPDIR="$outDir/$id.$hla_allele.$plen/tmp"
	mkdir -p $TMPDIR
	chmod 1777 $TMPDIR

	python /Share/BP/zhenglt/01.bin/netMHC/netMHC-3.4/netMHC.py -a "$hla_allele" -l $plen $pfile > "$outDir/netMHC.$id.$hla_allele.$plen.out"

	echo end at: `date` hla: $hla_allele len: $plen

	#rfile="/WPS1/zhenglt/work/MHC/TCGA/INPUT/TCGA.9.bed"
	#source /WPS1/zhenglt/work/MHC/TCGA/func.lib.sh
	#processNetMHCOut \\
	#	"$outDir/$id.9/$id.9.ori.out" \\
	#	"$outDir/$id.9/$id.9.mut.out" \\
	#	$rfile \\
	#	"$outDir/$id.9/$id.9.netMHC.out"
	#rm -r $outDir/$id.9/tmp
	#cut -f 12,13 "$outDir/$id.9/$id.9.netMHC.out.filtered.txt" | sort | uniq -c | awk -v OFS="\t" '{print \$2,\$3,\$1}' > "$outDir/$id.9/$id.9.netMHC.out.filtered.stat"
	#tmpfile="$outDir/$id.9/$id.9.netMHC.out.filtered.tmp"
	#listMerge.pl -i=0,1 -j=0,1 -a -b 0 \$tmpfile "/WPS1/zhenglt/work/MHC/TCGA/INPUT/TCGA.9.stat" "$outDir/$id.9/$id.9.netMHC.out.filtered.stat"
	#awk -F"\\t" -v OFS="\\t" 'BEGIN{print "CancerType\\tSample\\tNumOfMissense\\tNumOfBinder"} {print \$0}' \$tmpfile > "$outDir/$id.9/$id.9.netMHC.out.filtered.stat.txt"
	#rm "$outDir/$id.9/$id.9.netMHC.out.filtered.stat" \$tmpfile
	#/WPS1/zhenglt/work/MHC/TCGA/bin/binder.count.R "$outDir/$id.9/$id.9.netMHC.out.filtered.stat.txt" "$outDir/$id.9/$id.9.netMHC.out.filtered.stat"

	#NetMHCOutInd2Pop \\
	#	"$outDir/$id.9/$id.9.netMHC.out.filtered.txt" \\
	#	"/WPS1/zhenglt/work/MHC/TCGA/INPUT/merge.maf.used.out" \\
	#	"$outDir/$id.9/$id.9.netMHC.out.filtered.pop.txt"
	#NetMHCOutInd2Pop \\
	#	"$outDir/$id.9/$id.9.netMHC.out.txt" \\
	#	"/WPS1/zhenglt/work/MHC/TCGA/INPUT/merge.maf.used.out" \\
	#	"$outDir/$id.9/$id.9.netMHC.out.pop.txt"
	#cutIC50 "$outDir/$id.9/$id.9.netMHC.out.pop.txt" "$outDir/$id.9/$id.9.netMHC.out.pop.IC50"

}
