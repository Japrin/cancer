#!/bin/bash

sampleID="SAMPLE"
t1="File1"
t2="File2"
optL=0
optD=0
while getopts s:a:b:l:d: opt
do
	case $opt in 
	s)
		sampleID="$OPTARG"
		;;
	a)
		t1="$OPTARG"
		;;
	b)
		t2="$OPTARG"
		;;
	l)
		optL="$OPTARG"
		;;
	d)
		optD="$OPTARG"
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-s sampleID] [-a file1's title] [-b file2' title] [-l cnv size's cutoff] <file1 bed> <file2 bed/gff> <out prefix>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-s sampleID] [-a file1's title] [-b file2' title] [-l cnv size's cutoff] <file1 bed> <file2 bed/gff> <out prefix>"
	exit 1
fi

file1=$1
file2=$2
fileo=$3.overlap.bed
osummary=$3.overlap.summary

oDir=`dirname $3`
cmpfile1=$oDir/`basename $file1`.col3

sed 's/^chr//i' $file1 | awk -F"\t" -v OFS="\t" '{if($1=="23") $1="X"; if($1=="24") $1="Y"; if($3-$2>='$optL')print $1,$2,$3}' > $cmpfile1

if [[ "$file2" =~ ".gff" ]];then
	cmpfile2=$oDir/`basename $file2`
	echo "Format: gff"
	intersectBed -a $cmpfile1 -b $cmpfile2 -wo -f 0.6 -r | awk -F"\t" -v OFS="\t"  '{print $0,"'$sampleID'"}' > $fileo
	bpOverlap=`awk 'BEGIN{sum=0}{sum+=$(NF-1)}END{print sum}' $fileo`
	bpfile1=`awk 'BEGIN{sum=0}{sum+=$3-$2}END{print sum}' $cmpfile1`
	bpfile2=`awk 'BEGIN{sum=0}{sum+=$5-$4}END{print sum}' $cmpfile2`
	nOverlapFile1=`cut -f 1-3 $fileo | sort | uniq | wc -l`
	nOverlapFile2=`cut -f 4,7-8 $fileo | sort | uniq | wc -l`
elif [[ "$file2" =~ ".bed" ]];then
	cmpfile2=$oDir/`basename $file2`.col3
	awk -F"\t" -v OFS="\t" '{if($3-$2>='$optL') print $1,$2,$3}' $file2 > $cmpfile2
	echo "Format: bed"
	intersectBed -a $cmpfile1 -b $cmpfile2 -wo -f 0.5 -r | awk -F"\t" -v OFS="\t"  '{print $0,"'$sampleID'"}' > $fileo
	bpOverlap=`awk 'BEGIN{sum=0}{sum+=$(NF-1)}END{print sum}' $fileo`
	bpfile1=`awk 'BEGIN{sum=0}{sum+=$3-$2}END{print sum}' $cmpfile1`
	bpfile2=`awk 'BEGIN{sum=0}{sum+=$3-$2}END{print sum}' $cmpfile2`
	nOverlapFile1=`cut -f 1-3 $fileo | sort | uniq | wc -l`
	nOverlapFile2=`cut -f 4-6 $fileo | sort | uniq | wc -l`
fi

nTotalFile1=`cat $cmpfile1 | wc -l`
nTotalFile2=`cat $cmpfile2 | wc -l`
percentFile1=`echo "scale=4;100*$bpOverlap/$bpfile1" | bc`
percentFile2=`echo "scale=4;100*$bpOverlap/$bpfile2" | bc`
percentFile1ByNumber=`echo "scale=4;100*$nOverlapFile1/$nTotalFile1" | bc`
percentFile2ByNumber=`echo "scale=4;100*$nOverlapFile2/$nTotalFile2" | bc`

(
printf "Sample\ttotal_number_${t1}\ttotal_number_${t2}\toverlap_number_${t1}\toverlap_number_${t2}\tbase_number_${t1}\tbase_number_${t2}\toverlap_base_number\tsensitivity_by_number\tPPV_by_number\tsensitivity_by_base\tPPV_by_base\tcnv_size\tused_reads\n"
printf "$sampleID\t$nTotalFile1\t$nTotalFile2\t$nOverlapFile1\t$nOverlapFile2\t$bpfile1\t$bpfile2\t$bpOverlap\t$percentFile1ByNumber\t$percentFile2ByNumber\t$percentFile1\t$percentFile2\t$optL\t$optD\n"
)> $osummary
