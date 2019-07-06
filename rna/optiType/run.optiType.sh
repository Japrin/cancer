#!/bin/bash

optT="rna"
while getopts t: opt
do
	case $opt in 
	t)
		optT=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
        echo "Usage: $0 [-t 'rna'(default) or 'dna'] <outDir> <fq1> <fq2> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
    echo "Usage: $0 [-t 'rna'(default) or 'dna'] <outDir> <fq1> <fq2> <sampleID>"
	exit 1
fi

echo begin at: `date`

# export cbc
##export PATH=$PATH:/Share/BP/zhangfan/00.bin/Cbc-2.9.3/build/bin
export LD_LIBRARY_PATH=/WPSnew/zhenglt/01.bin/repseq/Cbc-2.9.3/build/lib:$LD_LIBRARY_PATH
export PATH=$PATH:/WPSnew/zhenglt/01.bin/repseq/Cbc-2.9.3/build/bin
### require package: pyomo
###export PYTHONPATH=/Share/BP/zhangfan/00.bin/python2.7/lib/python2.7/site-packages:$PYTHONPATH

##ot_dir="/Share/BP/zhangfan/05.pipeline/OptiType-fan"
ot_dir="/WPSnew/zhenglt/01.bin/repseq/OptiType-fan"
shDir=`dirname $(readlink -f $0)`

outDir=$1
fq1=$2
fq2=$3
sampleID=$4

mkdir -p $outDir
if [ "$optT" = "rna" ];then
    echo python $ot_dir/bin/optitype_pipeline_fan.py -o $outDir -c $shDir/config.rna.ini -i $fq1 $fq2
    python2 $ot_dir/bin/optitype_pipeline_fan.py -o $outDir -c $shDir/config.rna.ini -i $fq1 $fq2
else
    echo python $ot_dir/bin/optitype_pipeline_fan.py -o $outDir -c $shDir/config.dna.ini -i $fq1 $fq2
    python2 $ot_dir/bin/optitype_pipeline_fan.py -o $outDir -c $shDir/config.dna.ini -i $fq1 $fq2
fi

echo end at: `date`
