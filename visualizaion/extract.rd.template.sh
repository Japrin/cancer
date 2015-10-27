#!/bin/bash

iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"
binSize=100

while getopts c:b: opt
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
	b)
		binSize="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-b binSize, default 100] <sampleID> <infile> <outfile>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-b binSize, default 100] <sampleID> <infile> <outfile>"
	exit 1
fi

SHDIR=$(dirname $0)
export PATH=$SHDIR:$PATH

source $iniFile

sampleID=$1
infile=$2
ofile=$3
outDir=$(dirname $ofile)
#ofile=${infile/.root/.RD.gz}
#root_script="/PROJ/GR/share/medinfo.02pipeline/cancer/visualizaion/cnvnator_root2txt.C"
root_script="$SHDIR/cnvnator_root2txt.C"

echo begin at: `date` PID: $$
cd $outDir
#/PROJ/GR/share/Software/root/bin/root -b -q -l $root_script'("'$infile'",'$binSize')' \
#	| awk '/^chr/||/^#/' \
#	| bgzip -f -c \
#	> $ofile
#tabix -s 1 -b 2 -e 3 -f $ofile

if [ -f "/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/chr.24.length" ]
then
	while read chr length
	do
		$SHDIR/plot.cnvnator.region.sh $ofile chr$chr 1 $length $outDir/$sampleID.cnvnator.$chr 0
	done</PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/chr.24.length
fi

echo end at: `date`
