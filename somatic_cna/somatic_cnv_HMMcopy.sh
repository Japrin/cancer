#!/bin/bash

iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"
binSize=1000
afile="/PROJ/GR/share/medinfo.00database/genome/human/hg19/HMMcopy/hg19.gc.w1000.wig"
bfile="/PROJ/GR/share/medinfo.00database/genome/human/hg19/HMMcopy/hg19.map.w1000.wig"
while getopts c:w:a:b: opt
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
	w)
		binSize=$OPTARG
		;;
	a)
		if [ -f $OPTARG ]
		then
			afile="$OPTARG"
		else
			echo "WARNING: invalid file ($OPTARG), default will be used"
		fi
		;;
	b)
		if [ -f $OPTARG ]
		then
			bfile="$OPTARG"
		else
			echo "WARNING: invalid file ($OPTARG), default will be used"
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-w binSize, default 1000] [-a gc file] [-b map file] <sampleID> <normalBam> <tumorBam> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-w binSize, default 1000] [-a gc file] [-b map file] <sampleID> <normalBam> <tumorBam> <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile

sampleID=$1
normalBam=$2
tumorBam=$3
outDir=$4

mkdir -p $outDir
cd $outDir

if [ ! -f "$outDir/$sampleID.tumor_reads.wig" ];then
	/PROJ/GR/share/Software/medinfo/01bin/HMMcopy/HMMcopy/mybuild/bin/readCounter -w $binSize $normalBam > $outDir/$sampleID.normal_reads.wig
	/PROJ/GR/share/Software/medinfo/01bin/HMMcopy/HMMcopy/mybuild/bin/readCounter -w $binSize $tumorBam > $outDir/$sampleID.tumor_reads.wig
fi

/PROJ/GR/share/medinfo.02pipeline/cancer/somatic_cna/HMMcopy.core.R $outDir/$sampleID.normal_reads.wig $outDir/$sampleID.tumor_reads.wig $afile $bfile $outDir/$sampleID.HMMcopy $sampleID

### filter N
sed -e '1,1d' -e 's/^chr//' $outDir/$sampleID.HMMcopy.segment.txt \
	| awk -F"\t" -v OFS="\t" '$1~/^([0-9]|X)+/ && $4!=3{$2=$2-1;$3=$3+0;cn=$4-1;if(cn>2) svtype="gain"; else svtype="loss";print $1,$2,$3,$3-$2,svtype,"NA","NA","NA",cn,"NA"  }' \
	| sort -k 1,1 -k 2g,2 -k 3g,3 \
	> $outDir/$sampleID.HMMcopy.segment.CNV.bed

windowBed -a $outDir/$sampleID.HMMcopy.segment.CNV.bed -b /PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/all.NBlock.larger1000bp.bed -w 50000 -v > $outDir/$sampleID.HMMcopy.segment.CNV.final.bed

### pre-annotation (get .gff)
awk -F"\t" -v OFS="\t" 'BEGIN{id=0}{id++;print $1,"HMMcopy",$5,$2+1,$3,".",".",".","CopyNumber="$9";Size="$4";SVID="id";SVType="$5; print $1,"HMMcopy",$5,$2+1,$2+1,".",".",".","CopyNumber="$9";Size="$4";SVID="id";SVType=breakpoint"; print $1,"HMMcopy",$5,$3,$3,".",".",".","CopyNumber="$9";Size="$4";SVID="id";SVType=breakpoint";}' $outDir/$sampleID.HMMcopy.segment.CNV.final.bed > $outDir/$sampleID.HMMcopy.segment.CNV.final.gff

### annotation
/PROJ/GR/share/medinfo.02pipeline/cancer/all/var_annotation.sh $outDir/$sampleID.HMMcopy.segment.CNV.final.gff $sampleID

### plot
awk -F"\t" -v OFS="\t" 'BEGIN{print "Chromosome\tStart\tRatio\tMedianRatio\tCopyNumber"} NR>1{sub("chr","",$2);$3++;print $2,$3,2^($9),"NA",$10-1}' $outDir/$sampleID.HMMcopy.bin.txt > $outDir/$sampleID.HMMcopy.bin.forPlot.txt
/PROJ/GR/share/medinfo.02pipeline/cancer/somatic_cna/plot_CNA.byChr.R $outDir/$sampleID.HMMcopy.bin.forPlot.txt 2

echo end at: `date`
