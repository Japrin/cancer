#!/bin/bash

echo "*** somatic SNV comparison ***"

TR=""
optTR=""
iniFile="/WPS/GR/zhengliangtao/02pipeline/novo.med.seq/cancer/parameter/init_human.sh"

while getopts c:r: opt
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
	r)	
		if [ -f $OPTARG ]
		then
			TR="$OPTARG"
			optTR="-l $TR"
		else
			echo "WARNING: invalid target file ($OPTARG), no target will be used"
		fi
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] <outDir> <vcf1> <vcf2> <ID1> <ID2>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 5 ]
then 
	echo "Usage: $0 [-c iniFile] <outDir> <vcf1> <vcf2> <ID1> <ID2>"
	exit 1
fi

echo begin at: `date`

source $iniFile

outDir=$1
vcfFile1=$2
vcfFile2=$3
ID1=$4
ID2=$5

mkdir -p $outDir

venn3.pl -i 0,1 -n1 $ID1  -n2 $ID2 -title "Comparison of SNV called" $outDir $vcfFile1 $vcfFile2 > $outDir/${ID1}_${ID2}_venn3.log

echo end at: `date`
echo "*** Finished somatic SNV comparison ***"
