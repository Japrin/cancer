#!/bin/bash

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
opt_T="2"

while getopts c:t: opt
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
		opt_T=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-t type ,1 OR 2(default)] <infile> <id>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile] [-t type ,1 OR 2(default)] <infile> <id>"
	exit 1
fi

echo begin at: `date`

source $iniFile

infile=$1
id=$2

/Share/BP/zhenglt/02.pipeline/cancer/rna/run.mbased.R -t $opt_T $infile
awk -F"\t" -v OFS="\t" 'NR==1{print "sampleID",$0} NR>1{print "'$id'",$0}' $infile.ASE.${opt_T}s.Gene.filter > $infile.ASE.${opt_T}s.Gene.filter.txt
awk -F"\t" -v OFS="\t" 'NR==1{print "sampleID",$0} NR>1{print "'$id'",$0}' $infile.ASE.${opt_T}s.Gene.filter.isoform > $infile.ASE.${opt_T}s.Gene.filter.isoform.txt

echo end at: `date`
