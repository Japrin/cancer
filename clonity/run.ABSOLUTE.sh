#!/bin/bash

segData=$1
snvData=$2
outDir=$3
sampleID=$4

shDir=`dirname $0`
source "$shDir/../parameter/init_human.sh"

echo begin at: `date`
mkdir -p $outDir
cd $outDir


$shDir/run.ABSOLUTE.R $segData $snvData $outDir $sampleID

echo end at: `date`
