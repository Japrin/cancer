#!/bin/bash

function readcount
{
	sampleID=$1
	varfile=$2
	bamfile=$3
	outfile=$4
	
	. /usr/share/Modules/init/bash
	export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":$MODULEPATH
	module load novomedSeq/1.0
	module load samtools/0.1.19
	REF="/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/hg19.fa"
	
	echo begin at: `date`
	mkdir -p `dirname $outfile`	
	#awk -F"\t" -v OFS="\t" '{gsub("^chr","",$1);chr=$1;pos=$2; print "chr"chr,pos }' $varfile | sort -k 1,1 -k 2g,2 > $outfile.tmp0
	awk -F"\t" -v OFS="\t" '{gsub("^chr","",$1);chr=$1;beg=$2-120;end=$2+120; print "chr"chr,beg,end}' $varfile | sort -k 1,1 -k 2g,2 > $outfile.tmp1
	samtools mpileup -q 10 -Q 15 -f $REF -l $outfile.tmp1 $bamfile > $outfile.tmp2
	awk '$4>0' $outfile.tmp2 > $outfile.tmp3
	VarScan.jar readcounts $outfile.tmp3 --variants-file $outfile.tmp0 --output-file $outfile.tmp4 --min-coverage 1 --min-base-qual 15
	mv $outfile.tmp4 $outfile.RC
	rm $outfile.tmp0 $outfile.tmp1 $outfile.tmp2 $outfile.tmp3
	echo end at: `date`
}

readcount $1 $2 $3 $4
