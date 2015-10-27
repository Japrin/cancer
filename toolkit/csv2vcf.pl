#!/usr/bin/env perl
#============================================================================
# Name        		: filterAnnovarCSV.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Sat Jul 23 16:40:47 2011
# Last Modified By	: 
# Last Modified On	: Sat Jul 23 16:40:47 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl csv2vcf.pl [option] <infile>
	--type		GATK,samtools,varScan,with.format,other etc [defautl GATK]
	--sample	samples,comma seperated [required]
	--header	vcf header file [default build-in header]
	-h		display this help and exit
	
	Note: You must set --type be "other" if the file is neither from GATK nor samtools.

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use lib("/PUBLIC/software/HEALTH/04lib/perl/lib/perl5/site_perl/5.18.2");
use Text::CSV;


my ($in,$out);
my ($opt_h,$opt_type,$opt_samples,$opt_header);
GetOptions("h"	=>\$opt_h,
			"type=s"=>\$opt_type,
			"sample=s"=>\$opt_samples,
			"header=s"=>\$opt_header,
			);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;
$opt_type ||="GATK";
$opt_samples ||="sample";

my $csv = Text::CSV->new ( { binary => 1,eol => "\n" } ) or die "Cannot use CSV: ".Text::CSV->error_diag ();
if($infile=~/\.gz$/) { open $in,"gzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; } 
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
$_=<$in>;
chomp $_;
my @h=split /,/,$_;
for(my $i=0;$i<6;$i++) { pop @h; }
my @ss=split /,/,$opt_samples;
my %tagType=();
outputHeader($opt_type,$opt_header);
if($opt_type=~/GATK/i  || $opt_type=~/samtools/i || $opt_type=~/varScan/i  || $opt_type=~/with\.format/)
{
	printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",join("\t",@ss);
}else
{
	printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}
while ( my $row = $csv->getline($in) ) 
{
	#chr1	59374	.	A	G	255.07	HARD_TO_VALIDATE	AC=2;AF=1.00;AN=2;DP=20;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=0.9867;MQ=16.98;MQ0=8;QD=12.75;SB=-92.68	GT:AD:DP:GQ:PL	1/1:0,20:20:35.90:288,36,0
	#my @ii;
	my %ii;
	for(my $i=0;$i<@h;$i++)
	{
		if(!defined($row->[$i]) || $row->[$i] eq "") { next; }
		if($h[$i]=~/dbSNP/i) 
		{
			next; 
		}
		if($row->[$i]=~/;Name=/) { $row->[$i]=~s/;Name=/(/; $row->[$i].=")";}
		$row->[$i] =~ s/;/,/g;
		$ii{"$h[$i]=$row->[$i];"}=$i;
	}
	my @ii=sort {$ii{$a}<=>$ii{$b}} keys %ii;
	my $ii=sprintf "%s",join("",@ii);
	my $rsid=$row->[8]?$row->[8]:".";

	if($opt_type=~/GATK/i || $opt_type=~/samtools/i || $opt_type=~/varScan/i || $opt_type=~/with\.format/i)
	{
		my @_ss=();
		for(my $i=0;$i<@ss;$i++) { push @_ss,pop @$row; }
		my ($chr,$pos,$id,$ref,$alt,$qua,$filter,$info,$format)=@$row[-9 .. -1];
		#Func	Gene	ExonicFunc	AAChange	Conserved	SegDup	1000G_CEU	1000G_YRI	1000G_JPTCHB	dbSNP132	SIFT	PolyPhen2	LJB_PhyloP	LJB_MutationTaster	LJB_LRT	miRNA	tfbs	ensGene
		if($info eq ".") { $info=""; }
		printf "$chr\t$pos\t$rsid\t$ref\t$alt\t$qua\t$filter\t$ii$info\t$format\t%s\n",join("\t",(reverse @_ss));
	}else
	{
		while($row->[-1] eq "") { pop @$row; }
		my ($chr,$beg,$end,$ref,$alt,$h,$info)=@$row[-7 .. -1];
		print "$chr\t$beg\t$rsid\t$ref\t$alt\t.\t.\tzygosity=$h;$ii$info\n";
	}
	

}
$csv->eof or $csv->error_diag();
close $in;

############################################################################
sub usage
{
	die `pod2text $0`;
}
sub outputHeader
{
	my ($type,$_hFile)=@_;
	my $_header="";
	if($type=~/GATK/i)
	{
		if(defined($_hFile)) { $_header=`awk '/^##/' $_hFile`; }
		else
		{
			$_header=sprintf <<"Here";
##fileformat=VCFv4.1
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=StandardFilter,Description="QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=ABHet,Number=1,Type=Float,Description="Allele Balance for hets (ref/(ref+alt))">
##INFO=<ID=ABHom,Number=1,Type=Float,Description="Allele Balance for homs (A/(A+O))">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=OND,Number=1,Type=Float,Description="Overall non-diploid ratio (alleles/(alleles+non-alleles))">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
##INFO=<ID=Func,Number=.,Type=String,Description="Genome functional element(exonic,splicing,intergenic,intronic,ncRNA,upstream,downstream,UTR3,UTR5)">
##INFO=<ID=Gene,Number=.,Type=String,Description="Gene name">
##INFO=<ID=ExonicFunc,Number=1,Type=String,Description="exonic variant functoin(frameshift insertion, frameshift deletion, frameshift block substitution, stopgain, stoploss, nonframeshift insertion, nonframeshift deletion, nonframeshift block substitution, missense, synonymous, unknown)">
##INFO=<ID=1000g2012apr_ALL,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=ESP6500si_ALL,Number=1,Type=Float,Description="variant frequency in the ESP6500 population">
##INFO=<ID=cosmic64,Number=.,Type=String,Description="COSMIC">
##INFO=<ID=AAChange,Number=.,Type=String,Description="Amino acid change">
##INFO=<ID=mutType,Number=1,Type=String,Description="het(erozygous) or hom(ozygous)">
##INFO=<ID=Conserved,Number=1,Type=String,Description="Conserved score of the location region, range from 0 to 1000.(LOD socre of the conserved region is also present)">
##INFO=<ID=SegDup,Number=1,Type=String,Description="when variant locate in segmental duplications, this score represent sequence identity between two duplicate genome segement">
##INFO=<ID=AVSIFT,Number=1,Type=Float,Description="SIFT score">
##INFO=<ID=LJB_SIFT,Number=1,Type=Float,Description="LJB SIFT score">
##INFO=<ID=LJB_SIFT_Pred,Number=1,Type=String,Description="LJB SIFT prediction">
##INFO=<ID=LJB_PolyPhen2,Number=1,Type=Float,Description="PolyPhen2 score">
##INFO=<ID=LJB_PolyPhen2_Pred,Number=1,Type=String,Description="LJB PolyPhen2 prediction">
##INFO=<ID=LJB_PhyloP,Number=1,Type=Float,Description="PhyloP score">
##INFO=<ID=LJB_PhyloP_Pred,Number=1,Type=String,Description="LJB PhyloP prediction">
##INFO=<ID=LJB_MutationTaster,Number=1,Type=Float,Description="MutationTaster score">
##INFO=<ID=LJB_MutationTaster_Pred,Number=1,Type=String,Description="LJB mutationTaster prediction">
##INFO=<ID=LJB_LRT,Number=1,Type=Float,Description="LRT score">
##INFO=<ID=LJB_LRT_Pred,Number=1,Type=String,Description="LJB LRT prediction">
##INFO=<ID=LJB_GERP++,Number=1,Type=Float,Description="LJB GERP++">
##contig=<ID=1,assembly=b37,length=249250621>
##contig=<ID=2,assembly=b37,length=243199373>
##contig=<ID=3,assembly=b37,length=198022430>
##contig=<ID=4,assembly=b37,length=191154276>
##contig=<ID=5,assembly=b37,length=180915260>
##contig=<ID=6,assembly=b37,length=171115067>
##contig=<ID=7,assembly=b37,length=159138663>
##contig=<ID=8,assembly=b37,length=146364022>
##contig=<ID=9,assembly=b37,length=141213431>
##contig=<ID=10,assembly=b37,length=135534747>
##contig=<ID=11,assembly=b37,length=135006516>
##contig=<ID=12,assembly=b37,length=133851895>
##contig=<ID=13,assembly=b37,length=115169878>
##contig=<ID=14,assembly=b37,length=107349540>
##contig=<ID=15,assembly=b37,length=102531392>
##contig=<ID=16,assembly=b37,length=90354753>
##contig=<ID=17,assembly=b37,length=81195210>
##contig=<ID=18,assembly=b37,length=78077248>
##contig=<ID=19,assembly=b37,length=59128983>
##contig=<ID=20,assembly=b37,length=63025520>
##contig=<ID=21,assembly=b37,length=48129895>
##contig=<ID=22,assembly=b37,length=51304566>
##contig=<ID=X,assembly=b37,length=155270560>
##contig=<ID=Y,assembly=b37,length=59373566>
##contig=<ID=MT,assembly=b37,length=16569>
##contig=<ID=GL000207.1,assembly=b37,length=4262>
##contig=<ID=GL000226.1,assembly=b37,length=15008>
##contig=<ID=GL000229.1,assembly=b37,length=19913>
##contig=<ID=GL000231.1,assembly=b37,length=27386>
##contig=<ID=GL000210.1,assembly=b37,length=27682>
##contig=<ID=GL000239.1,assembly=b37,length=33824>
##contig=<ID=GL000235.1,assembly=b37,length=34474>
##contig=<ID=GL000201.1,assembly=b37,length=36148>
##contig=<ID=GL000247.1,assembly=b37,length=36422>
##contig=<ID=GL000245.1,assembly=b37,length=36651>
##contig=<ID=GL000197.1,assembly=b37,length=37175>
##contig=<ID=GL000203.1,assembly=b37,length=37498>
##contig=<ID=GL000246.1,assembly=b37,length=38154>
##contig=<ID=GL000249.1,assembly=b37,length=38502>
##contig=<ID=GL000196.1,assembly=b37,length=38914>
##contig=<ID=GL000248.1,assembly=b37,length=39786>
##contig=<ID=GL000244.1,assembly=b37,length=39929>
##contig=<ID=GL000238.1,assembly=b37,length=39939>
##contig=<ID=GL000202.1,assembly=b37,length=40103>
##contig=<ID=GL000234.1,assembly=b37,length=40531>
##contig=<ID=GL000232.1,assembly=b37,length=40652>
##contig=<ID=GL000206.1,assembly=b37,length=41001>
##contig=<ID=GL000240.1,assembly=b37,length=41933>
##contig=<ID=GL000236.1,assembly=b37,length=41934>
##contig=<ID=GL000241.1,assembly=b37,length=42152>
##contig=<ID=GL000243.1,assembly=b37,length=43341>
##contig=<ID=GL000242.1,assembly=b37,length=43523>
##contig=<ID=GL000230.1,assembly=b37,length=43691>
##contig=<ID=GL000237.1,assembly=b37,length=45867>
##contig=<ID=GL000233.1,assembly=b37,length=45941>
##contig=<ID=GL000204.1,assembly=b37,length=81310>
##contig=<ID=GL000198.1,assembly=b37,length=90085>
##contig=<ID=GL000208.1,assembly=b37,length=92689>
##contig=<ID=GL000191.1,assembly=b37,length=106433>
##contig=<ID=GL000227.1,assembly=b37,length=128374>
##contig=<ID=GL000228.1,assembly=b37,length=129120>
##contig=<ID=GL000214.1,assembly=b37,length=137718>
##contig=<ID=GL000221.1,assembly=b37,length=155397>
##contig=<ID=GL000209.1,assembly=b37,length=159169>
##contig=<ID=GL000218.1,assembly=b37,length=161147>
##contig=<ID=GL000220.1,assembly=b37,length=161802>
##contig=<ID=GL000213.1,assembly=b37,length=164239>
##contig=<ID=GL000211.1,assembly=b37,length=166566>
##contig=<ID=GL000199.1,assembly=b37,length=169874>
##contig=<ID=GL000217.1,assembly=b37,length=172149>
##contig=<ID=GL000216.1,assembly=b37,length=172294>
##contig=<ID=GL000215.1,assembly=b37,length=172545>
##contig=<ID=GL000205.1,assembly=b37,length=174588>
##contig=<ID=GL000219.1,assembly=b37,length=179198>
##contig=<ID=GL000224.1,assembly=b37,length=179693>
##contig=<ID=GL000223.1,assembly=b37,length=180455>
##contig=<ID=GL000195.1,assembly=b37,length=182896>
##contig=<ID=GL000212.1,assembly=b37,length=186858>
##contig=<ID=GL000222.1,assembly=b37,length=186861>
##contig=<ID=GL000200.1,assembly=b37,length=187035>
##contig=<ID=GL000193.1,assembly=b37,length=189789>
##contig=<ID=GL000194.1,assembly=b37,length=191469>
##contig=<ID=GL000225.1,assembly=b37,length=211173>
##contig=<ID=GL000192.1,assembly=b37,length=547496>
##contig=<ID=NC_007605,assembly=b37,length=171823>
##contig=<ID=hs37d5,assembly=b37,length=35477943>
Here
		}
		print $_header;
		my @_line=split /\n/,$_header;
		foreach (@_line)
		{
			if(/^##INFO/)
			{
				if(/ID=(.+?),.*Type=(.+?),/)
				{
					$tagType{$1}=$2;
				}else
				{
					die "bad header: $_\n";
				}
			}
		}
	}elsif($type=~/samtools/i)
	{
		if(defined($_hFile)) { $_header=`cat $_hFile`; }
		else
		{
			$_header=sprintf <<"Here";
##fileformat=VCFv4.1
##samtoolsVersion=0.1.18 (r982:295)
##INFO=<ID=Func,Number=.,Type=String,Description="Genome functional element(exonic,splicing,intergenic,intronic,ncRNA,upstream,downstream,UTR3,UTR5)">
##INFO=<ID=Gene,Number=.,Type=String,Description="Gene name">
##INFO=<ID=ExonicFunc,Number=1,Type=String,Description="exonic variant functoin(frameshift insertion, frameshift deletion, frameshift block substitution, stopgain, stoploss, nonframeshift insertion, nonframeshift deletion, nonframeshift block substitution, missense, synonymous, unknown)">
##INFO=<ID=1000g2012apr_ALL,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=ESP6500si_ALL,Number=1,Type=Float,Description="variant frequency in the ESP6500 population">
##INFO=<ID=AAChange,Number=.,Type=String,Description="Amino acid change">
##INFO=<ID=mutType,Number=1,Type=String,Description="het(erozygous) or hom(ozygous)">
##INFO=<ID=Conserved,Number=1,Type=String,Description="Conserved score of the location region, range from 0 to 1000.(LOD socre of the conserved region is also present)">
##INFO=<ID=SegDup,Number=1,Type=String,Description="when variant locate in segmental duplications, this score represent sequence identity between two duplicate genome segement">
##INFO=<ID=AVSIFT,Number=1,Type=Float,Description="SIFT score">
##INFO=<ID=LJB_SIFT,Number=1,Type=Float,Description="LJB SIFT score">
##INFO=<ID=LJB_SIFT_Pred,Number=1,Type=String,Description="LJB SIFT prediction">
##INFO=<ID=LJB_PolyPhen2,Number=1,Type=Float,Description="PolyPhen2 score">
##INFO=<ID=LJB_PolyPhen2_Pred,Number=1,Type=String,Description="LJB PolyPhen2 prediction">
##INFO=<ID=LJB_PhyloP,Number=1,Type=Float,Description="PhyloP score">
##INFO=<ID=LJB_PhyloP_Pred,Number=1,Type=String,Description="LJB PhyloP prediction">
##INFO=<ID=LJB_MutationTaster,Number=1,Type=Float,Description="MutationTaster score">
##INFO=<ID=LJB_MutationTaster_Pred,Number=1,Type=String,Description="LJB mutationTaster prediction">
##INFO=<ID=LJB_LRT,Number=1,Type=Float,Description="LRT score">
##INFO=<ID=LJB_LRT_Pred,Number=1,Type=String,Description="LJB LRT prediction">
##INFO=<ID=LJB_GERP++,Number=1,Type=Float,Description="LJB GERP++">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##FORMAT=<ID=PL,Number=.,Type=Integer,Description="List of Phred-scaled genotype likelihoods, number of values is (#ALT+1)*(#ALT+2)/2">
Here
		}
		print $_header;
		my @_line=split /\n/,$_header;
		foreach (@_line)
		{
			if(/^##INFO/)
			{
				if(/ID=(.+?),.*Type=(.+?),/)
				{
					$tagType{$1}=$2;
				}else
				{
					die "bad header: $_\n";
				}
			}
		}
	}elsif($type eq "other" || $type=~/with.format/)
	{
		if(defined($_hFile)) { $_header=`cat $_hFile`; }
		else
		{
			$_header=sprintf <<Here;
##fileformat=VCFv4.1
##INFO=<ID=Func,Number=.,Type=String,Description="Genome functional element(exonic,splicing,intergenic,intronic,ncRNA,upstream,downstream,UTR3,UTR5)">
##INFO=<ID=Gene,Number=.,Type=String,Description="Gene name">
##INFO=<ID=ExonicFunc,Number=1,Type=String,Description="exonic variant functoin(frameshift insertion, frameshift deletion, frameshift block substitution, stopgain, stoploss, nonframeshift insertion, nonframeshift deletion, nonframeshift block substitution, missense, synonymous, unknown)">
##INFO=<ID=1000g2012apr_ALL,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=ESP6500si_ALL,Number=1,Type=Float,Description="variant frequency in the ESP6500 population">
##INFO=<ID=AAChange,Number=.,Type=String,Description="Amino acid change">
##INFO=<ID=mutType,Number=1,Type=String,Description="het(erozygous) or hom(ozygous)">
##INFO=<ID=Conserved,Number=1,Type=String,Description="Conserved score of the location region, range from 0 to 1000.(LOD socre of the conserved region is also present)">
##INFO=<ID=SegDup,Number=1,Type=String,Description="when variant locate in segmental duplications, this score represent sequence identity between two duplicate genome segement">
##INFO=<ID=AVSIFT,Number=1,Type=Float,Description="SIFT score">
##INFO=<ID=LJB_SIFT,Number=1,Type=Float,Description="LJB SIFT score">
##INFO=<ID=LJB_SIFT_Pred,Number=1,Type=String,Description="LJB SIFT prediction">
##INFO=<ID=LJB_PolyPhen2,Number=1,Type=Float,Description="PolyPhen2 score">
##INFO=<ID=LJB_PolyPhen2_Pred,Number=1,Type=String,Description="LJB PolyPhen2 prediction">
##INFO=<ID=LJB_PhyloP,Number=1,Type=Float,Description="PhyloP score">
##INFO=<ID=LJB_PhyloP_Pred,Number=1,Type=String,Description="LJB PhyloP prediction">
##INFO=<ID=LJB_MutationTaster,Number=1,Type=Float,Description="MutationTaster score">
##INFO=<ID=LJB_MutationTaster_Pred,Number=1,Type=String,Description="LJB mutationTaster prediction">
##INFO=<ID=LJB_LRT,Number=1,Type=Float,Description="LRT score">
##INFO=<ID=LJB_LRT_Pred,Number=1,Type=String,Description="LJB LRT prediction">
##INFO=<ID=LJB_GERP++,Number=1,Type=Float,Description="LJB GERP++">
Here
		}
		print $_header;
	}elsif($type eq "varScan")
	{
		if(defined($_hFile)) { $_header=`cat $_hFile`; }
		else
		{
			$_header=sprintf <<Here;
##fileformat=VCFv4.1
##source=VarScan2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth of quality bases">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">
##INFO=<ID=SS,Number=1,Type=String,Description="Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)">
##INFO=<ID=SSC,Number=1,Type=String,Description="Somatic score in Phred scale (0-255) derived from somatic p-value">
##INFO=<ID=GPV,Number=1,Type=Float,Description="Fisher's Exact Test P-value of tumor+normal versus no variant for Germline calls">
##INFO=<ID=SPV,Number=1,Type=Float,Description="Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls">
##FILTER=<ID=str10,Description="Less than 10%% or more than 90%% of variant supporting reads on one strand">
##FILTER=<ID=indelError,Description="Likely artifact due to indel reads at this position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
##FORMAT=<ID=DP4,Number=.,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev">
##INFO=<ID=Func,Number=.,Type=String,Description="Genome functional element(exonic,splicing,intergenic,intronic,ncRNA,upstream,downstream,UTR3,UTR5)">
##INFO=<ID=Gene,Number=.,Type=String,Description="Gene name">
##INFO=<ID=ExonicFunc,Number=1,Type=String,Description="exonic variant functoin(frameshift insertion, frameshift deletion, frameshift block substitution, stopgain, stoploss, nonframeshift insertion, nonframeshift deletion, nonframeshift block substitution, missense, synonymous, unknown)">
##INFO=<ID=1000g2012apr_ALL,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=ESP6500si_ALL,Number=1,Type=Float,Description="variant frequency in the ESP6500 population">
##INFO=<ID=AAChange,Number=.,Type=String,Description="Amino acid change">
##INFO=<ID=Conserved,Number=1,Type=String,Description="Conserved score of the location region, range from 0 to 1000.(LOD socre of the conserved region is also present)">
##INFO=<ID=SegDup,Number=1,Type=String,Description="when variant locate in segmental duplications, this score represent sequence identity between two duplicate genome segement">
##INFO=<ID=AVSIFT,Number=1,Type=Float,Description="SIFT score">
##INFO=<ID=LJB_SIFT,Number=1,Type=Float,Description="LJB SIFT score">
##INFO=<ID=LJB_SIFT_Pred,Number=1,Type=String,Description="LJB SIFT prediction">
##INFO=<ID=LJB_PolyPhen2,Number=1,Type=Float,Description="PolyPhen2 score">
##INFO=<ID=LJB_PolyPhen2_Pred,Number=1,Type=String,Description="LJB PolyPhen2 prediction">
##INFO=<ID=LJB_PhyloP,Number=1,Type=Float,Description="PhyloP score">
##INFO=<ID=LJB_PhyloP_Pred,Number=1,Type=String,Description="LJB PhyloP prediction">
##INFO=<ID=LJB_MutationTaster,Number=1,Type=Float,Description="MutationTaster score">
##INFO=<ID=LJB_MutationTaster_Pred,Number=1,Type=String,Description="LJB mutationTaster prediction">
##INFO=<ID=LJB_LRT,Number=1,Type=Float,Description="LRT score">
##INFO=<ID=LJB_LRT_Pred,Number=1,Type=String,Description="LJB LRT prediction">
##INFO=<ID=LJB_GERP++,Number=1,Type=Float,Description="LJB GERP++">
Here
		}
		print $_header;
	}
}

