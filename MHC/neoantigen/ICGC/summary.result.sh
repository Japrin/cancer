#!/bin/bash

listFie=$1
outPrefix=$2

function summarizeOneSample
{
    annovarResFile=$1
    IEDBResFile=$2
    cancerType=$3
    sampleID=$4
    nMut=`cut -f 1-5 $annovarResFile | sort | uniq | wc -l`
    nEpi=`gzip -cd $IEDBResFile | awk -F"\t" -v OFS="\t" 'NR>1 && $28=="PASS" {print $19,$20,$21,$22,$23}'| sort | uniq  | wc -l`
    printf "$cancerType\t$sampleID\t$nMut\t$nEpi\n"
}

(
printf "CancerType\tSampleID\tMutBurden\tEpi\n"
while read cancerType sampleID IEDBResFile expFile annovarResFile
do
    summarizeOneSample $annovarResFile $IEDBResFile $cancerType $sampleID 
done<$listFie
) > $outPrefix.summary.tab

