#!/bin/bash -eu

iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"
optF=10000

while getopts c:f: opt
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
	f)
		optF="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-f flank,default 10000] <RDFile> <PEMFile> <libConfigFile> <chr> <begin> <end> <ID> <outDir>"
		printf "\n*Note*:\nRDFile: *.cnvnator.RD.txt.gz\nPEMFile: *.bam\nchr: must contain chr in the begining*\nlibConfigFile: with a header line:\n"
		printf "library\tuppercutoff\tlowercutoff\n"
		printf "NHD0657\t373.57\t144.35\n"
		printf "NHD0673\t389.53\t171.94\n"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 8 ]
then 
	echo "Usage: $0 [-c iniFile] [-f flank,default 10000] <RDFile> <PEMFile> <libConfigFile> <chr> <begin> <end> <ID> <outDir>"
	printf "\n*Note*:\nRDFile: *.cnvnator.RD.txt.gz\nPEMFile: *.bam\nchr: must contain chr in the begining*\nlibConfigFile: with a header line:\n"
	printf "library\tuppercutoff\tlowercutoff\n"
	printf "NHD0657\t373.57\t144.35\n"
	printf "NHD0673\t389.53\t171.94\n"
	exit 1
fi

echo begin at: `date`

source $iniFile 

RDFile=$1
PEMFile=$2
libConfigFile=$3
chr=$4
begin=$5
end=$6
SAMPLE=$7
outDir=$8

echo ">>> BEGIN visualization sv ..."

mkdir -p $outDir
aa=`echo "$begin-$optF" |bc`
bb=`echo "$end+$optF" |bc`

## RD input from cnvnator
/PROJ/GR/share/medinfo.02pipeline/cancer/visualizaion/plot.cnvnator.region.sh $RDFile $chr $begin $end $outDir/$SAMPLE $optF
## PEM input from breakdancer
/PROJ/GR/share/medinfo.02pipeline/cancer/visualizaion/breakdancer.abnormalRead.pl \
	-r $chr:${aa}-${bb} \
	-d \
	$libConfigFile \
	$PEMFile \
	>  $outDir/$SAMPLE.abnormalRead.txt \
	2> $outDir/$SAMPLE.abnormalRead.err
printf "$chr\t$aa\t$bb\n" | sed 's/^chr//' > $outDir/$SAMPLE.sv.bed
aa=`echo "$begin-1000" |bc`
bb=`echo "$begin+1000" |bc`
printf "$chr\t$aa\t$bb\n" | sed 's/^chr//' >> $outDir/$SAMPLE.sv.bed
aa=`echo "$end-1000" |bc`
bb=`echo "$end+1000" |bc`
printf "$chr\t$aa\t$bb\n" | sed 's/^chr//' >> $outDir/$SAMPLE.sv.bed

## mappability
gem-mappability-retriever query \
	/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/gem/human_g1k_v37_decoy.mappability.100bp.out.hash \
	$outDir/$SAMPLE.sv.bed \
	> $outDir/$SAMPLE.sv.mappability
## visualizaion	
perl /PROJ/GR/share/medinfo.02pipeline/cancer/visualizaion/sv.rd.pem.visualization.pl \
	-i $outDir/$SAMPLE.sv.bed \
	-o $outDir/$SAMPLE.sv \
	-b $outDir/$SAMPLE.abnormalRead.txt \
	-c $outDir/$SAMPLE.for.cnvnator.plot \
	-m $outDir/$SAMPLE.sv.mappability \
	-p $begin,$end

echo ">>> END visualization sv"
echo end at: `date`
