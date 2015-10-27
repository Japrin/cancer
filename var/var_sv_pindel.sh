#!/bin/bash -eu

iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"
opt_breakDancerFile=""

while getopts c:b: opt
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
		opt_breakDancerFile="-b $OPTARG"	
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-b breakdancer file] <sampleID> <outDir> <inbam>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-b breakdancer file] <sampleID> <outDir> <inbam>"
	exit 1
fi


source $iniFile 

#echo "*** Calling SV using pindel ***"
sampleID=$1
outDir=$2
inBam=$3

mkdir -p $outDir

confirmedBreakDancerFile=$outDir/$sampleID.confirmedBD.by.pindel.txt
isize=500
bamConfig=$outDir/$sampleID.pindel.bamConfig

echo begin at: `date`
printf "%s\t%s\t%s\n" $inBam $isize $sampleID >$bamConfig
cd $outDir
pindel -f $refData -i $bamConfig -o $outDir/$sampleID -c ALL -T 4 -d 85 -a 2 -g -L $outDir/$sampleID.pindel.log
#pindel -f $refData -i $bamConfig -o $outDir/$sampleID $opt_breakDancerFile -c ALL -T 4 -d 85 -a 2 -g -L $outDir/$sampleID.pindel.log

echo end at `date`
