#!/bin/bash

echo "*** Transform to svg file ***"

if [ $# -lt 2 ]
then
    echo "Usage: $0 <infile> <outSVG>"
    exit 1
fi

infile=$1
outSVG=$2
tmp=.vecj20121121lo

perl -ane 'if(/\s+(.+?)(\[.+)/){printf "    \"$1\"$2\n"}elsif(/\s+(.+?)\s+?->\s+(.+);/){printf "    \"$1\"  ->  \"$2\";\n";}else {print}' $infile > $tmp
dot -Tsvg $tmp > $outSVG
rm $tmp
