#!/bin/bash


iniFile="/WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/parameter/init_human.sh"
_refData="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"

optT=8
optA=""
optB=""
## another option: --solexa-quals ((-5, 40))
## another option: --phred64-quals (+64)

while getopts c:t:f:ab opt
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
		optT=$OPTARG
		;;
	f)
		if [ -f $OPTARG ]
		then
			_refData=$OPTARG
		fi
		;;
        a)
                #optA="--phred64-quals"
                optA="--solexa1.3-quals"
                ;;
        b)  
                optB="--fusion-search"
                ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	        echo "Usage: $0 [-c iniFile] [-t threads, default 8] [-a quality mode, +33 or +64 (default +33)] [-b fusion or not(default not)]  <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 6 ]
then 
	echo "Usage: $0 [-c iniFile] [-t threads, default 8] [-a quality mode, +33 or +64 (default +33)] [-b fusion or not(default not)]  <outDir> <fq1> <fq2, '-' if none> <sampleID> <lib> <lane> "
	exit 1
fi


source $iniFile
refData=$_refData
REF=$_refData

#genome="/WPS/BP/huboqiang/data/reference/homo_sapien/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
genome="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/bowtie2/human_g1k_v37_decoy"
annGTF="/WPS/BP/zhenglt/00.database/ensemble/release75/homo_sapiens/gtf/Homo_sapiens.GRCh37.75.gtf"
transcriptome="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/bowtie2/human_g1k_v37_decoy.transcriptome/Homo_sapiens.GRCh37.75"
IGVGenome="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_b37_decoy5.igv.genome"

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
lib=$5
lane=$6

mkdir -p $outDir

echo begin at: `date`
echo "*** Aligning reads ***"

#cd $outDir
   # -G $annGTF \

tophat \
    -o $outDir \
    -p $optT \
    --transcriptome-index $transcriptome \
    --library-type fr-unstranded \
    --read-mismatches  2 \
    --read-gap-length  2 \
    --read-edit-dist   2 \
    --rg-id $sampleID.$lib.$lane \
    --rg-sample  $sampleID \
    --rg-library $lib \
    --rg-platform-unit $lane \
    --keep-tmp \
    --keep-fasta-order \
    $genome $fq1 $fq2 $optA $optB 

samtools index $outDir/accepted_hits.bam

/WPS/BP/huboqiang/bin/igvtools count $outDir/accepted_hits.bam $outDir/$sampleID.tdf,$outDir/$sampleID.wig $IGVGenome

cufflinks --no-update-check -p $optT -u -G $annGTF -o $outDir $outDir/accepted_hits.bam
 
cuffquant --no-update-check -p $optT -o $outDir $annGTF $outDir/accepted_hits.bam 

 

echo end at: `date`
echo "*** Finished aligning reads ***"
