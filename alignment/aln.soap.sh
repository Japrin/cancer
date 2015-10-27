#!/bin/bash

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
faIndex="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/soap/human_g1k_v37_decoy.fasta.index"
fa="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/soap/human_g1k_v37_decoy.fasta"
NBlockFile="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/all.NBlock.larger1000bp.bed"
optT=4

while getopts c:t:f: opt
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
	t)
		optT=$OPTARG
		;;
	f)
		if [ -f $OPTARG ]
		then
			faIndex=$OPTARG
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-f reference] [-t threads, default 4] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-c iniFile] [-f reference] [-t threads, default 4] <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi

echo begin at: `date`

source $iniFile
export PATH="/Share/BP/zhenglt/01.bin/alignment/soap2/soap2.21release":$PATH

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
lib=$5
lane=$6

mkdir -p $outDir

soap -D $faIndex -a $fq1 -b $fq2 -o $outDir/$sampleID.soap.pair.ali -2 $outDir/$sampleID.soap.unpair.ali -u $outDir/$sampleID.soap.unmap.fa -p 4 -m 0 -x 2000

soap.coverage -cvg \
	-refsingle $fa \
	-i $outDir/$sampleID.soap.pair.alii $outDir/$sampleID.soap.unpair.ali \
	-o $outDir/$sampleID.soap.coverage.out \
	-depthsingle $outDir/$sampleID.soap.coverage.depthsingle.out \
       	-addn $NBlockFile \
	-plot $outDir/$sampleID.soap.coverage.plot.out 0 1000 \
	-window $outDir/$sampleID.soap.coverage.window1000.out 1000 

#########-depth $outDir/$sampleID.soap.coverage.depthDir \
echo end at: `date`
