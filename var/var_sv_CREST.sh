#!/bin/bash

iniFile="/PROJ/HEALTH/share/health.02pipeline/cancer/parameter/init_human.sh"
myhost=`hostname`
myport=5000
optG="M"
optS=""
opts=""
optE=""

_refData=/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta
_refData2bit=/PUBLIC/database/HEALTH/genome/human/b37_gatk/human_g1k_v37_decoy.2bit

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
		echo "Usage: $0 [-c iniFile] [-g gender default M] [-p port default 5000] [-r ref (fa)] [-b ref (2bit)] [-s range chr1 etc] [-e only \"extractSClip\" step] <sampleID> <outDir> <bam1(tumor)> [bam2(normal)]"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-g gender default M] [-p port default 5000] [-r ref (fa)] [-b ref (2bit)] [-s range chr1 etc] [-e only \"extractSClip\" step] <sampleID> <outDir> <bam1(tumor)> [bam2(normal)]"
	exit 1
fi

echo begin at: `date`

source $iniFile
### replace redData in iniFile
refData=$_refData
refData2bit=$_refData2bit
### other path
export PATH=/PROJ/HEALTH/share/01bin/CAP3/CAP3:$PATH
export PATH=/PROJ/HEALTH/share/01bin/blat/blatSrc35/mybuild:$PATH
export PATH=/PUBLIC/software/public/VarCall/crest:$PATH
CRESTDir=/PUBLIC/software/public/VarCall/crest

sampleID=$1
outDir=$2
tumorBam=$3
normalBam=$4
mkdir -p $outDir
cd $outDir
echo begin at: `date`

tumorBName=`basename $tumorBam`
optT="-d $tumorBam"
if [ -f "$normalBam" ]
then
	normalBName=`basename $normalBam`
	optT="-d $tumorBam -g $normalBam"
fi

if [ ! -f $outDir/$tumorBName$opts.cover ]
then
	perl $CRESTDir/extractSClip.pl -o $outDir -i $tumorBam -ref_genome $refData $optS
	if [ -f "$normalBam" ]
	then
		perl $CRESTDir/extractSClip.pl -o $outDir -i $normalBam -ref_genome $refData $optS
	fi
fi

if [ "$optE" == "E" ]
then
	echo extractSClip done at: `date`
	exit 0
fi
### setup gfServer
echo set up gfServer at host: $myhost

gfServer start $myhost $myport $refData2bit -canStop >$outDir/server.log 2>$outDir/server.err &

printf "query at: %s" `date` >$outDir/server.log
gfServer status $myhost $myport >>$outDir/server.log 2>>$outDir/server.log
while [[ "$?" != "0" ]]
do
	sleep 20
	printf "query at %s: " `date` >>$outDir/server.log
	gfServer status $myhost $myport >>$outDir/server.log 2>>$outDir/server.log
done

### CREST 
perl $CRESTDir/CREST.pl \
	-f $outDir/$tumorBName$opts.cover \
	$optT $optS \
	-p $tumorBName$opts \
	-o $outDir \
	--ref_genome $refData \
	-t $refData2bit \
	-blatserver $myhost \
	-blatport $myport 

### turn down server
gfServer stop $myhost $myport

### visulization
perl $CRESTDir/bam2html.pl \
	$optT \
	-i $outDir/$tumorBName$opts.predSV.txt \
	--ref_genome $refData \
	-o $outDir/$tumorBName$opts.predSV.html

### annotation
perl $PIPELINE/var/var_sv_CREST.toGff.pl $outDir/$tumorBName$opts.predSV.txt > $outDir/$tumorBName$opts.predSV.gff
var_annotation.sh -m CREST $outDir/$tumorBName$opts.predSV.gff $sampleID

echo end at: `date`
