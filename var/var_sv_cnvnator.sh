#!/bin/bash -eu

iniFile="/PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh"

binSize=100

while getopts b:c: opt
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
	b)
		binSize="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-b binSize, default 100] <sampleID> <outDir> <inbam>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-b binSize, default 100] <sampleID> <outDir> <inbam>"
	exit 1
fi

echo begin at: `date`

source $iniFile 

#echo "*** Calling CNV using Read-Depth Analysis: $CNVNATOR ***"
sampleID=$1
outDir=$2
inBam=$3

CNVNATOR="cnvnator"


p=$outDir/$sampleID.cnvnator
o=$p.gff

mkdir -p $outDir

echo "*******************************"
echo $LD_LIBRARY_PATH
echo "*******************************"

if [ ! -f $p.txt ]
then
	if [ ! -f $p.root ]
	then
		echo ">> extracting reads"
		$CNVNATOR -root $p.root -tree $inBam -unique
	fi
	echo ">> Generating histogram"
	$CNVNATOR -root $p.root -his $binSize -d `dirname $REF`/byChr
	echo ">> Calculating statistics"
	$CNVNATOR -root $p.root -stat $binSize
	echo ">> RD signal partitioning"
	$CNVNATOR -root $p.root -partition $binSize
	echo ">> CNV calling"
	$CNVNATOR -root $p.root -call $binSize > $p.raw
	cat $p.raw | grep -v "\(WARN\)\|\(==\)" > $p.txt
fi

echo ">> Converting output to GFF"

NFile=/PUBLIC/database/HEALTH/genome/human/b37_gatk/all.NBlock.larger20bp.bed

minsize=0
awkopt='{feat=$4; if ($4=="deletion") feat="Deletion"; else if ($4=="duplication") feat="Duplication"; if ($5>='$minsize') print $1"\tCNVnator\t"feat"\t"$2"\t"$3"\t.\t.\t.\tSize="$5";RD="$6";q0="$7};'

ifile=$p.txt
ofile=$o

perl -F"\t" -ane 'chomp @F;$F[1]=~s/^chr//;$F[1]=~s/[:-]/\t/g;print join("\t",$F[1],$F[0],$F[2],$F[3],$F[8])."\n";' $ifile \
	| awk "$awkopt" | sort -k1,1 -k4n -k5n -u \
	| awk -v OFS="\t" '{print $1,$4-1,$5,$0}' \
	| intersectBed -a stdin -b $NFile -f 0.5 -r -v \
	| perl -F"\t" -ane 'chomp @F;shift @F;shift @F;shift @F;printf "%s\n",join("\t",@F);' \
	| awk '$1~/^[0-9XY]+$/' \
	| perl -F"\t" -ane '/Size=(.+?);RD=(.+?);q0=(.+)/;if(/^#/|| ($1>=1000 && ($2<0.6 || $2>1.4 ) && $3<0.5 && $F[0]!~/[Y]/) ){print "$_" }' \
	| perl -ane 'BEGIN{$SVID=0} $SVID++; chomp; @F=split /\t/; print "$_;SVID=$SVID;SVType=$F[2]\n";' \
	> $ofile

var_annotation.sh -c $iniFile $ofile $sampleID

echo end at: `date`
echo "*** Finished CNV Calling using Read-Depth Analysis"
