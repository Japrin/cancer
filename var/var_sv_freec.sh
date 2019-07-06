#!/bin/bash -eu

iniFile="`dirname $0`/../parameter/init_human.sh"

export freecDir="/WPSnew/zhenglt/01.bin/cancer/cna/freec/FREEC-11.5"
optG="/WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/gem/human_g1k_v37_decoy_150.mappability"
optS="/WPSnew/zhenglt/01.bin/cancer/cna/freec/hg19_snp142.SingleDiNucl.1based.txt"
optF="/WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/byChr"
optL="/WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.chr.length.h24"
optE="/WPSnew/zhenglt/02.pipeline/cancer/var/freec.config.example"
optP=2
optM="FR"
while getopts c:g:s:f:l:e:p:m: opt
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
	m)
		optM=$OPTARG
		;;
	g)
		optG=$OPTARG
		;;
	s)
		optS=$OPTARG
		;;
	f)
		optF=$OPTARG
		;;
	l)
		optL=$OPTARG
		;;
	e)
		optE=$OPTARG
		;;
	p)
		optP=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-g mappability file] [-s SNPfile] [-m mateOrientation] [-e template config] [-f *.fa dir] [-l length file] [-p ploidy, default 2] <sampleID> <outDir> <bam1(tumor)> [bam2(normal)]"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-g mappability file] [-s SNPfile] [-m mateOrientation] [-e template config] [-f *.fa dir] [-l length file] [-p ploidy, default 2] <sampleID> <outDir> <bam1(tumor)> [bam2(normal)]"
	exit 1
fi

echo begin at: `date`

source $iniFile 

sampleID=$1
shift
outDir=$1
shift

p=$outDir/`basename $1`

bams=''
for i in $*
do
	bams="$bams $i"
done

mkdir -p $outDir

if [ ! -f ${p}_CNVs ]
then
	if [ ! -f $outDir/config_freec ]
	then
		perl $PIPELINE/cancer/var/freec.generate.config.pl \
			-o $outDir \
			-gemMappabilityFile $optG -SNPfile $optS -chrFiles $optF -chrLenFile $optL -e $optE -m $optM \
			$bams \
			> $outDir/config_freec
	fi
    $freecDir/src/freec -conf $outDir/config_freec > ${p}_log.txt 2> ${p}_err.txt
    ##cat $freecDir/scripts/makeGraph.R | R --slave --args $optP ${p}_ratio.txt
    
    if [ ! -f ${p}_BAF.txt ]
	then
		cat $freecDir/scripts/makeGraph.R | R --slave --args $optP ${p}_ratio.txt
	else
		cat $freecDir/scripts/makeGraph.R | R --slave --args $optP ${p}_ratio.txt ${p}_BAF.txt
	fi
fi

function plotByChr
{
    local iChr=$1
    local optP=$2
    local p=$3
    cat $freecDir/scripts/makeGraph_Chromosome.R | R --slave --args $iChr $optP ${p}_ratio.txt ${p}_BAF.txt
    ##echo "cat $freecDir/scripts/makeGraph_Chromosome.R | R --slave --args $iChr $optP ${p}_ratio.txt ${p}_BAF.txt"
}
export -f plotByChr

parallel -j 8 plotByChr {} ::: {1..22} X Y ::: $optP ::: $p

exit

nCol=`awk 'NR==1{print NF}' ${p}_CNVs`
cp ${p}_CNVs ${p}_CNVs.copy
mv ${p}_CNVs ${p}_CNVs.bak
if [ $nCol -gt 5 ]
then	
	awk '$4!=2 && /somatic/' ${p}_CNVs.bak | $PIPELINE/cancer/somatic_cna/freec.cnv.region.pl > ${p}_CNVs	## somatic 
	$PIPELINE/cancer/somatic_cna/freec.cnv.region.pl ${p}_CNVs.bak > ${p}_CNVs.profile
	awk -v OFS="\t" '$1~/^([0-9]|X|Y)/{print $1,$2,$3,$3-$2,$5,"NA","NA","NA",$4,$6}' ${p}_CNVs > ${p}_CNVs.bed
else
	awk '$4!=2 ' ${p}_CNVs.bak > ${p}_CNVs
	awk -v OFS="\t" '$1~/^([0-9]|X|Y)/{print $1,$2,$3,$3-$2,$5,"NA","NA","NA",$4,"NA"}' ${p}_CNVs > ${p}_CNVs.bed
fi
### filter N
###ln -s -f ${p}_CNVs.bed ${p}_CNVs.final.bed	
intersectBed -a ${p}_CNVs.bed -b /WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/all.NBlock.larger1000bp.bed -f 0.5 -v > ${p}_CNVs.final.bed 
###windowBed -a ${p}_CNVs.bed -b /PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/all.NBlock.larger1000bp.bed -w 50000 -v > ${p}_CNVs.final.bed

### pre-annotation (get .gff)
awk -F"\t" -v OFS="\t" 'BEGIN{id=0}{id++;print $1,"FREEC",$5,$2,$3,".",".",".","CopyNumber="$9";Size="$4";SVID="id";SVType="$5; print $1,"FREEC",$5,$2,$2,".",".",".","CopyNumber="$9";Size="$4";SVID="id";SVType=breakpoint"; print $1,"FREEC",$5,$3,$3,".",".",".","CopyNumber="$9";Size="$4";SVID="id";SVType=breakpoint";}' ${p}_CNVs.final.bed > ${p}_CNVs.final.gff

### summary
$PIPELINE/cancer/somatic_cna/freec.cnv.summary.pl -s $sampleID ${p}_CNVs > ${p}_CNVs.summary.txt

if [ -f ${p}_BAF.txt ]
then
	awk -F"\t" -v OFS="\t" 'NR>1{$2=$2-1"\t"$2;print}' ${p}_BAF.txt \
		| intersectBed -a stdin -b ${p}_CNVs.profile -loj \
		| cut -f 2 --complement \
		| $PIPELINE/cancer/somatic_cna/freec.loh.region.pl \
		> ${p}_LOH.bd.marker.bed
	cut -f 1-5 --complement ${p}_LOH.bd.marker.bed | uniq > ${p}_LOH.bd.cnv.bed
	
	$PIPELINE/cancer/somatic_cna/freec.loh.summary.pl -s $sampleID ${p}_LOH.bd.marker.bed >${p}_tmp1
	$PIPELINE/cancer/somatic_cna/freec.cnv.summary.pl -s $sampleID ${p}_LOH.bd.cnv.bed    >${p}_tmp2
	printf "sample\tCNV_somatic_status\tCN\tcount(bd.markder)\tlength(bd.marker)\tcount(bd.cnv)\tlength(bd.cnv)\n" > ${p}_LOH.summary.txt
	paste ${p}_tmp1 ${p}_tmp2 | cut -f 1-5,9,10 >> ${p}_LOH.summary.txt
	rm ${p}_tmp1 ${p}_tmp2
fi
### annotation
var_annotation.sh ${p}_CNVs.final.gff $sampleID
#### plot
#perl $PIPELINE/var/freec.for.plot.pl ${p}_ratio.txt > ${p}_ratio.forPlot.bin.txt 2> ${p}_ratio.forPlot.seg.txt


echo end at: `date`
