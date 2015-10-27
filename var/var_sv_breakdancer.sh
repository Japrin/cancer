#!/bin/bash -eu

echo "*** Calling SV using Read-Pair Mapping: ***"

iniFile="/PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh"

optT=""
optO=""
optG="M"
optN=6

while getopts c:g:to:n: opt
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
	g)
		optG="$OPTARG"
		;;
	t)	
		optT="-t"
		;;
	n)
		optN=$OPTARG
		;;
	o)
		optO="-o $OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG"
		echo "Usage: $0 [-c iniFile] [-g gender default M] [-n supporting reads per sample, default 6] [-t only inter-chromosomal translocation] [-o chr] <sampleID> <outDir> <bam1(tumor)> [bam2(normal) ...]"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-g gender default M] [-n supporting reads per sample, default 6] [-t only inter-chromosomal translocation] [-o chr] <sampleID> <outDir> <bam1(tumor)> [bam2(normal) ...]"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
shift
outDir=$1
shift

bams=''
for i in $*
do
	bams="$bams `cd \`dirname $i\`; pwd`/`basename $i`"
done

echo ">> Generating configuration file"

mkdir -p $outDir

p=$outDir/$sampleID.breakdancer
o=$outDir/$sampleID.breakdancer.gff

cd $outDir


if [ ! -f $p.txt ]
then
	perl $BREAKDANCER/bam2cfg.pl -g -h $bams > $p.cfg
	echo ">> Performing read-pair mapping"
	$BREAKDANCER/breakdancer-max $optO $optT -h -d $p.SV-supporting -g $p.bed $p.cfg > $p.txt
fi

#mkdir -p $outDir/tigra_dump
#/WPS/GR/zhengliangtao/01bin/tigra-sv/mybin/bin/tigra-sv -b $p.txt $bams \
#	-R $refData \
#	-r -d \
#	-I $outDir/tigra_dump \
#	> $outDir/$sampleID.breakdancer.tigra.txt

## filter and format
var_sv_breakdancer.filter.pl -g $optG -n $optN -a $p.txt > $p.flt.txt
var_sv_breakdancer.toGff.pl $p.flt.txt >  $p.flt.gff

## annotation
var_annotation.sh -m BreakDancer $p.flt.gff $sampleID

## stat
breakdancer.sv.stat.pl -s $sampleID $p.flt.gff.ann.xls > $p.flt.gff.ann.stat.txt

## by sv type
#awk 'NR==1||/Deletion/' $outDir/$sampleID.BreakDancer.ann.txt > $outDir/$sampleID.BreakDancer.Deletion.ann.txt
#awk 'NR==1||/Inversion/' $outDir/$sampleID.BreakDancer.ann.txt > $outDir/$sampleID.BreakDancer.Inversion.ann.txt
#awk 'NR==1||/Translocation/' $outDir/$sampleID.BreakDancer.ann.txt > $outDir/$sampleID.BreakDancer.Translocation.ann.txt
#awk 'NR==1||/Insertion/' $outDir/$sampleID.BreakDancer.ann.txt > $outDir/$sampleID.BreakDancer.Insertion.ann.txt

## compare with DGV (Deletion)
#head -1 $outDir/$sampleID.BreakDancer.Deletion.ann.txt > $outDir/$sampleID.BreakDancer.Deletion.novel.ann.txt
#intersectBed -b /PROJ/GR/share/medinfo.00database/DGV/GRCh37_hg19_variants_2013-07-23.deletion.bed -a $outDir/$sampleID.BreakDancer.Deletion.ann.txt -wa -r -f 0.5 -v >> $outDir/$sampleID.BreakDancer.Deletion.novel.ann.txt
#breakdancer.sv.stat.pl -s $sampleID $outDir/$sampleID.BreakDancer.Deletion.novel.ann.txt > $outDir/$sampleID.BreakDancer.Deletion.novel.ann.stat.txt

#total=`awk '!/^#/' $outDir/$sampleID.BreakDancer.Deletion.ann.txt | wc -l`
#novel=`awk '!/^#/' $outDir/$sampleID.BreakDancer.Deletion.novel.ann.txt | wc -l`
#printf "===========================================\n" >> $outDir/$sampleID.BreakDancer.Deletion.novel.ann.stat.txt
#printf "#sample\tvarType\ttotal\tnovel\trecall rate\n" >> $outDir/$sampleID.BreakDancer.Deletion.novel.ann.stat.txt
#printf "$sampleID\tDeletion\t$total\t$novel\t%4.2f\n" `echo "scale=6;1-$novel/$total"|bc` >> $outDir/$sampleID.BreakDancer.Deletion.novel.ann.stat.txt

echo "*** Finished Calling SV using Read-Pair Mapping ***"
echo end at: `date`
