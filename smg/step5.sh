#!/bin/bash

ofile=all.coverage.roi.txt
printf "sample\tCovered\tATs_Covered\tCpGs_Covered\tCGs_Covered\n" > $ofile
while read id
do
	tail -n 1 coverage/roi_covgs/$id.covg | awk -v OFS="\t" '{print "'$id'",$2,$3,$4,$5}'

done<S30.list >> $ofile

ofile=all.coverage.roi.mutationRate.txt
paste <(sort -k 1r,1 all.coverage.roi.txt) <(sort -k 1r,1 roi.mut.txt) | cut -f 6 --complement > $ofile
