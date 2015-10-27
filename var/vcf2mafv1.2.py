import re
import gzip
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--tumor",
                  action="store", type="string", dest="tvcf_file",help="the tumor vcf.gz file")
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
#title = ['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_Position','End_Position','Strand','Variant_Classification',
#         'Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode',
#         'Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2',
#         'Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status','Mutation_Status','Sequencing_Phase',
#         'Sequence_Source','Validation_Method','Sequencer','Tumor_Sample_UUID','Matched_Norm_Sample_UUID','t_ref_count','t_alt_count']

title = ['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_position','End_Position','Strand','Variant_Classification',
         'Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode',
         'Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2',
         'Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status','Mutation_Status','Sequencing_Phase',
         'Sequence_Source','Validation_Method','Sequencer','Tumor_Sample_UUID','Matched_Norm_Sample_UUID','t_ref_count','t_alt_count']

B.writelines('\t'.join(title)+'\n')
var_class = {'missense_SNV':'Missense_Mutation',
             'synonymous_SNV':'Silent',
             'stopgain_SNV':'Nonsense_Mutation',
             'stoploss_SNV':'Nonstop_Mutation',
             'unknown':'unknown',
             'exonic,splicing':'Splice_Site',
             'nonframeshift_deletion':'In_Frame_Del',
             'frameshift_deletion':'Frame_Shift_Del',
             'nonframeshift_insertion':'In_Frame_Ins',
             'frameshift_insertion':'Frame_Shift_Ins',
             'ncRNA_UTR5':'RNA',
             'ncRNA_UTR3':'RNA',
             'UTR5':'5\'UTR',
             'UTR3':'3\'UTR',
	     'UTR5,UTR3':'5\'UTR',
             'intronic':'Intron',
             'ncRNA_intronic':'RNA',
             'ncRNA_exonic':'RNA',
             'ncRNA_splicing':'RNA',
	     'intergenic':'Intergenic',
             'upstream':'5\'Flank',
             'downstream':'3\'Flank',
             'upstream,downstream':'5\'Flank',
             '':''}

    
A = safe_open(tvcf,'r')
p = re.compile('Func=(.*?);')
q = re.compile('Gene=(.*?);')
r = re.compile('ExonicFunc=(.*?);')

for n in A:
    if n.startswith('#'):continue
    lines = []
    l = n.strip().split()
    m = q.findall(n)
    if m:a = m[0]
    lines.append(a)
    if a in gene:lines.append(gene[a])
    else:lines.append('')
    lines +=['biopic','hg19']
	#chromosome, beg, end
    l[0]=re.sub(r'','',l[0])
    lines += [l[0],l[1],l[1]]
    #strand
    lines.append('+')
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
    #if 'insertion' in b: lines.append('INS')
    #elif 'deletion' in b: lines.append('DEL')
    #else:lines.append('SNP')
    ref=l[3]
    alt=l[4]
    if len(ref) < len(alt) : lines.append('INS')
    elif len(ref) > len(alt) : lines.append('DEL')
    else: lines.append('SNP')
    #lines.append('SNP')
    lines.append(l[3])
    if l[10].split(':')[0] == '0/1':lines += [l[3],l[4],l[2],'']
    else:lines += [l[4],l[4],l[2],'']
    #Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode
    lines += ['','']
    lines += [l[3],l[3]]
    x = 'Somatic'
    
    lines += ['','','','','Unknown','Untested',x,'',sourc,'none','Illumina HiSeq','','']
    lines += [l[10].split(':')[3],l[10].split(':')[4]]
    lines[15] = sampleID
    B.writelines('\t'.join(lines)+'\n')
A.close()
B.close()

