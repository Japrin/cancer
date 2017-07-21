#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optT=4
optS="human"
optM=""
optN=8

while getopts c:t:s:m: opt
do
	case $opt in 
	t)
		optT=$OPTARG
		;;
    s)
        optS=$OPTARG
        ;;
    m)
        optM=$OPTARG
        ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-t threads, default 4] [-s species, default human] [-m mode, default ''] <outDir> <fq1> <fq2, '-' if none> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-t threads, default 4] [-s species, default human] [-m mode, default ''] <outDir> <fq1> <fq2, '-' if none> <sampleID>"
	exit 1
fi

source $iniFile
module load gcc/4.9.2
module load kallisto/0.43.1
module load trimmomatic/0.33

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
TMPDIR="/Share/BP/zhenglt/tmp"
if [ "$optS" == "human" ];then
    indexFile="/DBS/DB_temp/zhangLab/ensemble/mybuild/kallisto/b20170519/rel88/Homo_sapiens.ensemble.release88.GRCh38.cdna.lincRNA.noCHR"
elif [ "$optS" == "mouse" ];then
    indexFile="/DBS/DB_temp/zhangLab/ensemble/mybuild/kallisto/b20170715/rel89/mouse/Mus_musculus.ensemble.release89.GRCm38.cdna.lincRNA.noCHR"
fi
optOO=""
if [ "$fq2" == "-" ];then
    fq2="--single"
    optOO="--fragment-length=200 --sd=60"
fi

echo begin at: `date`
mkdir -p $outDir

newFQ1=$fq1
newFQ2=$fq2

if [ "$optM" == "trim" ];then
	### filter low quality reads
    if [ "$fq2" == "--single" ];then
	    java -Xmx${optN}g -jar $TRIMMOMATIC_DIR/trimmomatic-0.33.jar SE \
            -threads 8 \
            $fq1 \
            $outDir/$sampleID.clean.fq.gz \
            ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-SE.fa:2:30:10 \
		    CROP:75 \
		    TRAILING:3 \
		    MAXINFO:50:0.25 \
		    MINLEN:36 \
		    TOPHRED33
        newFQ1=$outDir/$sampleID.clean.fq.gz
        newFQ2="--single"
    else
		    ##ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-PE-2.fa:2:30:10 \
	    java -Xmx${optN}g -jar $TRIMMOMATIC_DIR/trimmomatic-0.33.jar PE \
		    -threads 8 \
		    $fq1 \
		    $fq2 \
		    $outDir/$sampleID.clean.P.R1.fq.gz \
		    $outDir/$sampleID.clean.UP.R1.fq.gz \
		    $outDir/$sampleID.clean.P.R2.fq.gz \
		    $outDir/$sampleID.clean.UP.R2.fq.gz \
		    ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/NexteraPE-PE.fa:2:30:10 \
		    CROP:75 \
		    TRAILING:3 \
		    MAXINFO:50:0.25 \
		    MINLEN:36 \
		    TOPHRED33 
        newFQ1=$outDir/$sampleID.clean.P.R1.fq.gz
        newFQ2=$outDir/$sampleID.clean.P.R2.fq.gz
    fi
	echo "trimmomatic done. (" `date` ")"
fi	

echo kallisto quant --index=$indexFile -o $outDir --bootstrap-samples=100 --threads=$optT $newFQ1 $newFQ2 $optOO
kallisto quant --index=$indexFile -o $outDir --bootstrap-samples=100 --threads=$optT $newFQ1 $newFQ2 $optOO

echo end at: `date`
