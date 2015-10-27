#!/bin/env python

import sys
import os
import re

if len(sys.argv) <= 2:
	print 'usage: <output file> <gff file>'
	exit(1)
	
path = os.environ['ANNOVAR']
dPath = os.environ['HumanDB']
output = sys.argv[1]
input = sys.argv[2]

if input.endswith('.gff'):
	gff = input
	avinput = gff + '.avinput'
	avoutput = gff + '.avoutput'
	gff_file = open(gff, 'read')
	out_file = open(avinput, 'w')
	for line in gff_file:
		if re.match('^#',line) != None:
			continue
		gffCols = line.rstrip('\n').split('\t')
		out_file.write(gffCols[0]+'\t'+gffCols[3]+'\t'+gffCols[4]+'\t0\t0\t'+gffCols[2]+'\t'+gffCols[5]+'\t'+gffCols[8]+'\n')
	out_file.flush()
	out_file.close()
else:
	print >> sys.stderr, "Unknown input format: "+input
	exit(1)

# use user-define output name for now
avoutput=output

print 'Annotating variants with hg19 UCSC refGene...\n'
os.system('%s/annotate_variation.pl --geneanno --buildver hg19 -dbtype refGene --separate %s %s' %(path, avinput, dPath))
print 'Annotating variants with hg19 RMSK...\n'
os.system('%s/annotate_variation.pl -regionanno --buildver hg19 -dbtype gff3 -gff3dbfile hg19_rmsk.gff %s %s' %(path, avinput, dPath))
print 'Annotating variants with hg19 segment duplication...\n'
os.system('%s/annotate_variation.pl -regionanno --buildver hg19 -dbtype genomicSuperDups %s %s' %(path, avinput, dPath) )

exonic_file = avinput + '.exonic_variant_function'
function_file = avinput + '.variant_function'
rmsk_file = avinput + ".hg19_gff3"
segdup_file= avinput + ".hg19_genomicSuperDups"

def makeDict(filename, chrCol, startCol, endCol, valueCols, isGFF=False):
	file = open(filename, 'r')
	dic = {}
	for line in file:
		cols = line.split('\t')
		key=(cols[chrCol], cols[startCol], cols[endCol])
		if key not in dic:
			dic[key] = []
			for valueCol in valueCols:
				dic[key].append([])
		for i in range(len(valueCols)):
			value=cols[valueCols[i]]
			if isGFF:
				value=value.split(";")[-1].split("=")[-1]
			dic[key][i].append(value)
	return dic

function = makeDict(function_file, 2, 3, 4, [1, 0])
rmsk = makeDict(rmsk_file, 2, 3, 4, [1], True)
segdup = makeDict(segdup_file,2,3,4,[1])

AVINPUT = open(avinput, 'r')
AVOUTPUT = open(avoutput, 'w')
AVOUTPUT.write("#chr\tstart\tend\tgene_symbol\tregion\trmsk\tsegmental duplication\tother_info");

AVOUTPUT.write("\n")

def write(dic, key, nvalueCols=1):
	if key in dic: 
		for value in dic[key]:
			AVOUTPUT.write('\t'+";".join(value))
	else:		
		AVOUTPUT.write('\t.'*nvalueCols)

for line in AVINPUT:
	#if re.match('([0-9A-Za-z]+)\s+(\d+)', line):
		splitline = line.rstrip('\n').split('\t')
		key=(splitline[0], splitline[1],splitline[2])
		AVOUTPUT.write(splitline[0]+'\t'+splitline[1]+'\t'+splitline[2])
		write(function, key, 2)
		write(rmsk, key)
		write(segdup, key)
		
		AVOUTPUT.write('\t'+splitline[5]+'\t'+splitline[6]+'\t'+splitline[7])
		AVOUTPUT.write('\n')
	#else:
		#AVOUTPUT.write('Error\t'+line)

AVINPUT.flush();
AVOUTPUT.flush();

AVINPUT.close();
AVOUTPUT.close();

#`rm $avinput.*`
