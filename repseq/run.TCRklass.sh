#!/bin/bash

if [ $# -lt 4 ]
then 
	echo "Usage: $0 <outDir> <sampleID> <fq1> <fq2>"
	exit 1
fi

#outDir=testOut
#fq1=hs1a_1.fastq.gz
#fq2=hs1a_2.fastq.gz

outDir=$1
sampleID=$2
fq1=$3
fq2=$4

mkdir -p $outDir/TRA
mkdir -p $outDir/TRB

export PATH=/Share/BP/zhenglt/01.bin/RepSeq/TCRklass/htqc-0.90.8-Source/mybuild/bin:$PATH
export PERL5LIB=/Share/BP/zhenglt/04.lib/perl/lib/perl5/x86_64-linux-thread-multi:/Share/BP/zhenglt/04.lib/perl/lib/perl5:/Share/BP/zhenglt/04.lib/perl/lib64/perl5:/Share/BP/zhenglt/04.lib/perl/share/perl5

. /usr/share/Modules/init/bash
export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
module load TCRklass/0.6.0
module load boost/1.47

TCRHumanDBDir="/Share/BP/zhenglt/01.bin/RepSeq/TCRklass/benchmark_data/db/human"

echo begin at: `date`

echo .......... prepare input ..........
if [[ "$fq1" =~ ".gz" ]]; then
	gzip -cd $fq1 > $outDir/$sampleID.R1.fq
	gzip -cd $fq2 > $outDir/$sampleID.R2.fq
	fq1=$outDir/$sampleID.R1.fq
	fq2=$outDir/$sampleID.R2.fq
fi

echo .......... analyze TRA ..........
tcrklass.pl --in $fq1 $fq2 --out $outDir/TRA --species human \
	--v-nt $TCRHumanDBDir/av_n.db \
	--j-nt $TCRHumanDBDir/aj_n.db \
	--v-aa $TCRHumanDBDir/av_p.db \
	--j-aa $TCRHumanDBDir/aj_p.db

echo .......... analyze TRB ..........
tcrklass.pl --in $fq1 $fq2 --out $outDir/TRB --species human \
	--v-nt $TCRHumanDBDir/bv_n.db \
	--j-nt $TCRHumanDBDir/bj_n.db \
	--v-aa $TCRHumanDBDir/bv_p.db \
	--j-aa $TCRHumanDBDir/bj_p.db

if [[ "$fq1" =~ ".gz" ]]; then
	rm $outDir/$sampleID.R1.fq
	rm $outDir/$sampleID.R2.fq
fi
echo end at: `date`
