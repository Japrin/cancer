#!/bin/bash

iniFile="/PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh"
myhost=`hostname`
myport=5000
optG="M"
optS=""
opts=""
optE=""

_refData=/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/human_g1k_v37_decoy.fasta
_refData2bit=/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/human_g1k_v37_decoy.2bit

while getopts c:p:g:r:b:s:e opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid ini file ($OPTARG), default will be used"
		fi
		;;
	s)
		optS="-r $OPTARG"
		opts=".$OPTARG"
		;;
	e)
		optE="E"
		;;
	g)
		optG="$OPTARG"
		;;
	p)
		myport=$OPTARG	
		;;
	r)
		if [ -f $OPTARG ]
		then
			_refData="$OPTARG"
		else
			echo "WARNING: invalid ref(fa) file ($OPTARG), default will be used"
		fi
		;;
	b)
		if [ -f $OPTARG ]
		then
			_refData2bit="$OPTARG"
		else
			echo "WARNING: invalid ref(2bit) file ($OPTARG), default will be used"
		fi
		;;
	'?')
		echo "Usage: $0 [-c iniFile] [-g gender default M] [-r ref (fa)] [-b ref (2bit)] <sampleID> <outDir> <bam1(tumor)> [bam2(normal)]"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-g gender default M] [-r ref (fa)] [-b ref (2bit)] <sampleID> <outDir> <bam1(tumor)> [bam2(normal)]"
	exit 1
fi

echo begin at: `date`

source $iniFile
### replace redData in iniFile
refData=$_refData
refData2bit=$_refData2bit
### other path
DellyDir=/PROJ/HEALTH/share/01bin/delly/delly/mybin
PATH=$DellyDir:$PATH
export PATH=$PATH

sampleID=$1
outDir=$2
tumorBam=$3
normalBam=$4
mkdir -p $outDir
cd $outDir
echo begin at: `date`

delly -t DEL -x $DellyDir/human.hg19.excl.tsv -q 30 -o $outDir/$sampleID.delly.DEL.vcf -g $refData $tumorBam $normalBam
delly -t DUP -x $DellyDir/human.hg19.excl.tsv -q 30 -o $outDir/$sampleID.delly.DUP.vcf -g $refData $tumorBam $normalBam
delly -t INV -x $DellyDir/human.hg19.excl.tsv -q 30 -o $outDir/$sampleID.delly.INV.vcf -g $refData $tumorBam $normalBam
delly -t TRA -x $DellyDir/human.hg19.excl.tsv -q 30 -o $outDir/$sampleID.delly.TRA.vcf -g $refData $tumorBam $normalBam

### annotation
#perl $PIPELINE/var/var_sv_CREST.toGff.pl $outDir/$tumorBName$opts.predSV.txt > $outDir/$tumorBName$opts.predSV.gff
#var_annotation.sh -m CREST $outDir/$tumorBName$opts.predSV.gff $sampleID

echo end at: `date`
