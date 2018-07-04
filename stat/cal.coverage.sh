#!/bin/bash

TR=""
optTR=""
sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"
statBinDir="$sDir/../stat"
####iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
_refData=""

while getopts c:f:r: opt
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
			optTR="-L $TR"
		else
			echo "WARNING: invalid target file ($OPTARG), no target will be used"
		fi
		;;
	f)
		_refData=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r targetRegion] [-f refData] <sampleID> <inBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-r targetRegion] [-f refData] <sampleID> <inBam> <outDir>"
	exit 1
fi

source "$iniFile"

echo begin at: `date`

sampleID=$1
inBam=$2
outDir=$3

if [ -f "$_refData" ]
then
	refData=$_refData
	REF=$refData
fi

mkdir -p $outDir
java -Djava.io.tmpdir=$outDir -Xmx5g -jar $gatkJAR -T DepthOfCoverage \
	-mmq 0 \
	-mbq 0 \
	-omitBaseOutput \
	-R $REF \
	-I $inBam $optTR \
	-o $outDir/$sampleID \
	-et NO_ET \
	-K $gatkKey \
	-rf BadCigar \
	--interval_merging OVERLAPPING_ONLY \
	-ct 1 \
	-ct 4 \
	-ct 10 \
	-ct 15 \
	-ct 20 \
	-ct 30 \
	-ct 40 \
	-ct 50 \
	-ct 100

perl -i -F"\t" -ane '{ $NR++;if($NR==1){print; next} chomp; chomp @F; @a=split /[:-]/,$F[0];if(@a<3){push @a,$a[1];} $F[0]="$a[0]:$a[1]-$a[2]";$F[0]=~s/^chr//i;printf "%s\n",join("\t",@F); }' $outDir/$sampleID.sample_interval_summary 
#sed -i 's/^chr//' $outDir/$sampleID.sample_interval_summary
$statBinDir/raw.depth.gatk.wrapper.sh $outDir/$sampleID

$statBinDir/CovByChr.sh $outDir/$sampleID.sample_interval_summary $outDir/$sampleID.coverage.bychr.txt

perl -F"\t" -ane 'chomp @F;$F[0]=~s/^chr//;if($F[0] ne "Y"){$sum+=$F[2];$count+=$F[1];}if($F[0] eq "Y"){$YDep=$F[3]} END{ if($YDep<$sum/(3*$count)){printf "gender=F\n"; }else{ printf "gender=M\n";} }' $outDir/$sampleID.coverage.bychr.txt > $outDir/$sampleID.gender

echo end at: `date`
