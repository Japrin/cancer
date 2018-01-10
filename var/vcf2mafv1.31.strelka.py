#!/usr/bin/env python


import re
import gzip
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--tumor",
                  action="store", type="string", dest="tvcf_file",help="the tumor vcf.gz file OR '-'(from stdin)")
parser.add_option("-n", "--normal",
                  action="store", type="string", dest="nvcf_file",help="the normal vcf.gz file",default="")
parser.add_option("-s", "--seq_source",
                  action="store", type="string", dest="seq_source",help="WGS or other")
parser.add_option("-g", "--gene_data",
                  action="store", type="string", dest="gene_file",help="gene_name & Entrez_Gene_Id")

parser.add_option("-o", "--output",
                  action="store", type="string", dest="maf_file",help="the maf file")
parser.add_option("-a", "--sample",
                  action="store", type="string", dest="sampleID",help="sampleID")

options, args = parser.parse_args()

def safe_open(file_name,mode):
    try:
        if not file_name.endswith('.gz'):
            return open(file_name,mode)
        else:
            import gzip
            return gzip.open(file_name,mode)
    except IOError:
        print file_name + ' do not exist!'


tvcf = options.tvcf_file
nvcf = options.nvcf_file
g_file = options.gene_file
maf = options.maf_file
sourc = options.seq_source
sampleID = options.sampleID

A = safe_open(g_file,'r')
gene = {}
for n in A:
    l = n.strip().split()
    gene[l[0]] = l[1]
A.close()

B = safe_open(maf,'w')

title = ['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_position','End_Position','Strand','Variant_Classification',
         'Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode',
         'Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2',
         'Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status','Mutation_Status','Sequencing_Phase',
         'Sequence_Source','Validation_Method','Sequencer','Tumor_Sample_UUID','Matched_Norm_Sample_UUID','t_ref_count','t_alt_count']

B.writelines('\t'.join(title)+'\n')
var_class = {'missense_SNV':'Missense_Mutation',
             'synonymous_SNV':'Silent',
             'stopgain_SNV':'Nonsense_Mutation',
             'stopgain':'Nonsense_Mutation',
             'stoploss_SNV':'Nonstop_Mutation',
             'stoploss':'Nonstop_Mutation',
             'unknown':'unknown',
             'nonframeshift_deletion':'In_Frame_Del',
             'frameshift_deletion':'Frame_Shift_Del',
             'nonframeshift_insertion':'In_Frame_Ins',
             'frameshift_insertion':'Frame_Shift_Ins',
			 'frameshift_substitution':'Frame_Shift_Ins',
             'ncRNA_UTR5':'RNA',
             'ncRNA_UTR3':'RNA',
             'UTR5':'5\'UTR',
             'UTR3':'3\'UTR',
			 'UTR5,UTR3':'5\'UTR',
             'intronic':'Intron',
             'ncRNA_intronic':'RNA',
             'ncRNA_exonic':'RNA',
			 'ncRNA_splicing':'RNA',
			 'exonic,splicing':'Splice_Site',
			 'ncRNA_exonic,splicing':'Splice_Site',
             'intergenic':'IGR',
             'upstream':'5\'Flank',
             'downstream':'3\'Flank',
             'upstream,downstream':'5\'Flank',
             '':''}


if tvcf=="-":
    A = sys.stdin
else:
    A = safe_open(tvcf,'r')
p = re.compile('Func=(.*?);')
q = re.compile('Gene=(.*?);')
r = re.compile('ExonicFunc=(.*?);')
#aaChangeP = re.compile('AAChange=ARID1A:uc001bmt.1:exon1:c.C770T:p.A257V,ARID1A:uc001bmu.1:exon1:c.C770T:p.A257V,ARID1A:uc001bmv.1:exon1:c.C770T:p.A257V;')

for n in A:
    if n.startswith('#'):continue
	###INFO=<ID=QSS,Number=1,Type=Integer,Description="Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal">
	##INFO=<ID=TQSS,Number=1,Type=Integer,Description="Data tier used to compute QSS">
	##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
	##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
	##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description="Data tier used to compute QSS_NT">
	##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
	##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
	##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
	##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
	##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
	##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
	##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
	##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
	##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
	##FILTER=<ID=BCNoise,Description="Fraction of basecalls filtered at this site in either sample is at or above 0.4">
	##FILTER=<ID=SpanDel,Description="Fraction of reads crossing site with spanning deletions in either sample exceeeds 0.75">
	##FILTER=<ID=QSS_ref,Description="Normal sample is not homozygous ref or ssnv Q-score < 15, ie calls with NT!=ref or QSS_NT < 15">
	#1       17084464        .       G       A       .       PASS    Func=exonic;Gene=MST1L;ExonicFunc=missense_SNV;AAChange=MST1L:uc001azp.5:exon6:c.C434T:p.S145F,MST1L:uc010ock.3:exon12:c.C1634T:p.S545F;phastConsElements46way=(Score=638,Name=lod=522);genomicSuperDups=(Score=0.987113,Name=chr1:16972126);cytoband=1p36.13;avsift=0;ljb2_pp2hdiv=1.0;ljb2_pp2hvar=1.0;NT=ref;QSS=102;QSS_NT=102;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1  DP:FDP:SDP:SUBDP:AU:CU:GU:TU    348:0:0:0:0,0:0,0:348,398:0,0   208:0:0:0:11,12:0,0:197,216:0,0
	#16      46552384        .       C       CA      .       PASS    IC=9;IHP=9;NT=ref;QSI=48;QSI_NT=48;RC=8;RU=A;SGT=ref->het;SOMATIC;TQSI=2;TQSI_NT=2      DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50	47:47:41,45:0,0:4,4:47.77:0.00:0.00     28:28:14,16:9,9:3,3:32.39:0.00:0.00
    lines = []
    l = n.strip().split()
    m = q.findall(n)
    if m:a = m[0]
    l[4]=l[4].split(',')[0]
    lines.append(a)
    if a in gene:lines.append(gene[a])
    else:lines.append('')
    lines +=['BIOPIC','hg19']
	#chromosome, beg, end
    l[0]=re.sub(r'chr','',l[0])
    lines += [l[0],l[1],l[1]]
    #strand
    lines.append('+')
    ## func
    m = p.findall(n)
    if m:
        a = m[0]
    if a == 'exonic':
        m = r.findall(n)
        if m:
            b = m[0]
            lines.append(var_class[b])
    elif a == 'splicing':lines.append('Splice_Site')
    else:lines.append(var_class[a])

    vType = ""
    if len(l[3]) < len(l[4]):
        lines.append('INS')
        vType = "INDEL"
    elif len(l[3]) > len(l[4]):
        lines.append('DEL')
        vType = "INDEL"
    else:
        lines.append('SNP')
        vType = "SNP"

    lines.append(l[3])
    ### for SNP, het or hom
    ###if l[10].split(':')[0] == '0/1':lines += [l[3],l[4],l[2],'']
    ###else:lines += [l[4],l[4],l[2],'']
    ### for somatic mutation, always "het"
    lines += [l[3],l[4],l[2],'']

    #Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode
    lines += ['','']
    lines += [l[3],l[3]]
    x = 'Somatic'

    lines += ['','','','','Unknown','Untested',x,'',sourc,'none','Illumina HiSeq','','']
    ref_dp = 0
    alt_dp = 0

	#1       17084464        .       G       A       .       PASS    Func=exonic;Gene=MST1L;ExonicFunc=missense_SNV;AAChange=MST1L:uc001azp.5:exon6:c.C434T:p.S145F,MST1L:uc010ock.3:exon12:c.C1634T:p.S545F;phastConsElements46way=(Score=638,Name=lod=522);genomicSuperDups=(Score=0.987113,Name=chr1:16972126);cytoband=1p36.13;avsift=0;ljb2_pp2hdiv=1.0;ljb2_pp2hvar=1.0;NT=ref;QSS=102;QSS_NT=102;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1  DP:FDP:SDP:SUBDP:AU:CU:GU:TU    348:0:0:0:0,0:0,0:348,398:0,0   208:0:0:0:11,12:0,0:197,216:0,0
	#16      46552384        .       C       CA      .       PASS    IC=9;IHP=9;NT=ref;QSI=48;QSI_NT=48;RC=8;RU=A;SGT=ref->het;SOMATIC;TQSI=2;TQSI_NT=2      DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50	47:47:41,45:0,0:4,4:47.77:0.00:0.00     28:28:14,16:9,9:3,3:32.39:0.00:0.00
    if vType == 'SNP':
        if l[3] == "A":
            ref_dp=(l[10].split(':')[4]).split(',')[1]
        elif l[3] == "C":
            ref_dp=(l[10].split(':')[5]).split(',')[1]
        elif l[3] == "G":
            ref_dp=(l[10].split(':')[6]).split(',')[1]
        elif l[3] == "T":
            ref_dp=(l[10].split(':')[7]).split(',')[1]

        if l[4] == "A":
            alt_dp=(l[10].split(':')[4]).split(',')[1]
        elif l[4] == "C":
            alt_dp=(l[10].split(':')[5]).split(',')[1]
        elif l[4] == "G":
            alt_dp=(l[10].split(':')[6]).split(',')[1]
        elif l[4] == "T":
            alt_dp=(l[10].split(':')[7]).split(',')[1]
    else:
        alt_dp=(l[10].split(':')[3]).split(',')[1]
        ref_dp=str(int(l[10].split(':')[0])-int(alt_dp))
    lines += [ref_dp,alt_dp]
    lines[15] = sampleID
    B.writelines('\t'.join(lines)+'\n')
A.close()
B.close()
