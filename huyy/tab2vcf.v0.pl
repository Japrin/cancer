#!/usr/bin/perl -w

=pod

=head1 NAME

			tab2vcf.pl

=head1 SYNOPSIS

perl tab2vcf.pl <infile> -id [ID] -format [FORMAT]

=head1 DESCRIPTION

reformat tab delimitated file(after reformat_annovar.pl) to vcf file

=head1 OPTIONS

=over

=item -id

specify the ID of the tab file

=item -format

FORMAT

=back

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($id,$format,$h,$m); 
GetOptions(
	"id=s"=>\$id,
	"format=s"=>\$format,
	"h|help"=>\$h,
	"m|man"=>\$m
);

if ($h) {pod2usage(-verbose => 1);}
if ($m) {pod2usage(-verbose  => 2);}
if (!defined $h and !defined $m and !defined $id and !defined $format)  { pod2usage(-verbose=>0);}

my $filename = shift;
open In, $filename or die "$!";

my @id = split ",", $id;
my $id_length = @id;

while (<DATA>){
	print $_;
}

print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
foreach (@id) {
	print "\t$_";
}
print "\n";


my (@headers,%headers);
while (<In>){
	chomp;
	my @record = split "\t";
	if ($record[-1] eq "Shared") {
		pop @record;
	}
	if (/^CHROM/) {
		foreach (0 .. $#record) {
			$record[$_] =~ s/\..+$//;
			$headers{$record[$_]} = $_;
			push @headers, $record[$_];
		}
		
	}
	else {
		my $dbsnp = $record[$headers{"snp137"}];
		foreach (0..6) {
			if (($_ == 2 ) and ($dbsnp)) {
				print "$dbsnp\t";
				next;
			}
			print "$record[$_]\t";
		}
		my $infoend = $#headers - 1 - $id_length;
		my @information;
		foreach (7..$infoend){
			next if ($record[$_] eq ".");
			next if ($headers[$_] eq "INFO");
			$record[$_] =~ s/ /_/g;
			$record[$_] =~ s/;/,/g;
			if ($record[$_] =~ /=/){
				$record[$_] = "("."$record[$_]".")";
			}
			push @information, "$headers[$_]=$record[$_]";
		}
		my $information = join ";", @information;
		print $information;
		my $info = $record[$headers{"INFO"}];
		if ($info) {
			print ";$info";
		}
		$format = $record[$headers{"FORMAT"}];
		print "\t$format";
		foreach (@id) {
			print "\t$record[$headers{$_}]";
		}
		print "\n";
	}
}

__DATA__
##fileformat=VCFv4.1
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=StandardFilter,Description="QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0">
##FILTER=<ID=PASS,Description="High quality">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality non-reference bases">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=Gene,Number=.,Type=String,Description="Gene name">
##INFO=<ID=Func,Number=.,Type=String,Description="Genome functional element(exonic,splicing,intergenic,intronic,ncRNA,upstream,downstream,UTR3,UTR5)">
##INFO=<ID=ExonicFunc,Number=1,Type=String,Description="exonic variant functoin(frameshift insertion, frameshift deletion, frameshift block substitution, stopgain, stoploss, nonframeshift insertion, nonframeshift deletion, nonframeshift block substitution, missense, synonymous, unknown)">
##INFO=<ID=1000g2012apr_all,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=esp6500si_all,Number=1,Type=Float,Description="variant frequency in the ESP6500 population">
##INFO=<ID=AAChange,Number=.,Type=String,Description="Amino acid change">
##INFO=<ID=cpgIslandExt,Number=.,Type=String,Description="CpG Island">
##INFO=<ID=evofold,Number=.,Type=String,Description="Conserve RNA prediction">
##INFO=<ID=wgRna,Number=.,Type=String,Description="snoRNA and miRNA annotation">
##INFO=<ID=targetScanS,Number=.,Type=String,Description="TargetScan generated miRNA target site predictions">
##INFO=<ID=phastConsElements46way,Number=.,Type=String,Description="conserved elements produced by the phastCons program based on a whole-genome alignment of vertebrates">
##INFO=<ID=genomicSuperDups,Number=.,Type=String,Description="Segmental duplications in genome">
##INFO=<ID=tfbsConsSites,Number=.,Type=String,Description="transcription factor binding sites conserved in the human/mouse/rat alignment, based on transfac Matrix Database (v7.0)">
##INFO=<ID=gwasCatalog,Number=.,Type=String,Description="Published GWAS results on diverse human diseases">
##INFO=<ID=cytoband,Number=.,Type=String,Description="the approximate location of bands seen on Giemsa-stained chromosomes">
##INFO=<ID=snp137NonFlagged,Number=.,Type=String,Description="dbSNP version 137 with ANNOVAR index files, after removing those flagged SNPs (SNPs < 1% minor allele frequency (MAF) (or unknown), mapping only once to reference assembly, flagged in dbSnp as 'clinically associated')">
##INFO=<ID=snp137,Number=.,Type=String,Description="dbSNP version 137">
##INFO=<ID=targetScanS,Number=.,Type=String,Description="miRNA target prediction by TargetScan">
##INFO=<ID=avsift,Number=1,Type=Float,Description="SIFT score">
##INFO=<ID=ljb2_pp2hdiv,Number=1,Type=Float,Description="whole-exome PolyPhen scores built on HumanDiv database (for complex phenotypes)">
##INFO=<ID=ljb2_pp2hvar,Number=1,Type=String,Description="whole-exome PolyPhen version 2 scores built on HumanVar database (for Mendelian phenotypes)">
##INFO=<ID=ljb_pp2,Number=1,Type=String,Description="whole-exome PolyPhen scores">
##INFO=<ID=gerp++elem,Number=1,Type=Float,Description="conserved genomic regions by GERP++">
##INFO=<ID=cosmic65,Number=.,Type=String,Description="description of cosmic">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=IS,Number=2,Type=Float,Description="Maximum number of reads supporting an indel and fraction of indel reads">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##INFO=<ID=QBD,Number=1,Type=Float,Description="Quality by Depth: QUAL/#reads">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Read Position Bias">
##INFO=<ID=MDV,Number=1,Type=Integer,Description="Maximum number of high-quality nonRef reads in samples">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias (v2) for filtering splice-site artefacts in RNA-seq data. Note: this version may be broken.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
##contig=<ID=GL000207.1,length=4262,assembly=b37>
##contig=<ID=GL000226.1,length=15008,assembly=b37>
##contig=<ID=GL000229.1,length=19913,assembly=b37>
##contig=<ID=GL000231.1,length=27386,assembly=b37>
##contig=<ID=GL000210.1,length=27682,assembly=b37>
##contig=<ID=GL000239.1,length=33824,assembly=b37>
##contig=<ID=GL000235.1,length=34474,assembly=b37>
##contig=<ID=GL000201.1,length=36148,assembly=b37>
##contig=<ID=GL000247.1,length=36422,assembly=b37>
##contig=<ID=GL000245.1,length=36651,assembly=b37>
##contig=<ID=GL000197.1,length=37175,assembly=b37>
##contig=<ID=GL000203.1,length=37498,assembly=b37>
##contig=<ID=GL000246.1,length=38154,assembly=b37>
##contig=<ID=GL000249.1,length=38502,assembly=b37>
##contig=<ID=GL000196.1,length=38914,assembly=b37>
##contig=<ID=GL000248.1,length=39786,assembly=b37>
##contig=<ID=GL000244.1,length=39929,assembly=b37>
##contig=<ID=GL000238.1,length=39939,assembly=b37>
##contig=<ID=GL000202.1,length=40103,assembly=b37>
##contig=<ID=GL000234.1,length=40531,assembly=b37>
##contig=<ID=GL000232.1,length=40652,assembly=b37>
##contig=<ID=GL000206.1,length=41001,assembly=b37>
##contig=<ID=GL000240.1,length=41933,assembly=b37>
##contig=<ID=GL000236.1,length=41934,assembly=b37>
##contig=<ID=GL000241.1,length=42152,assembly=b37>
##contig=<ID=GL000243.1,length=43341,assembly=b37>
##contig=<ID=GL000242.1,length=43523,assembly=b37>
##contig=<ID=GL000230.1,length=43691,assembly=b37>
##contig=<ID=GL000237.1,length=45867,assembly=b37>
##contig=<ID=GL000233.1,length=45941,assembly=b37>
##contig=<ID=GL000204.1,length=81310,assembly=b37>
##contig=<ID=GL000198.1,length=90085,assembly=b37>
##contig=<ID=GL000208.1,length=92689,assembly=b37>
##contig=<ID=GL000191.1,length=106433,assembly=b37>
##contig=<ID=GL000227.1,length=128374,assembly=b37>
##contig=<ID=GL000228.1,length=129120,assembly=b37>
##contig=<ID=GL000214.1,length=137718,assembly=b37>
##contig=<ID=GL000221.1,length=155397,assembly=b37>
##contig=<ID=GL000209.1,length=159169,assembly=b37>
##contig=<ID=GL000218.1,length=161147,assembly=b37>
##contig=<ID=GL000220.1,length=161802,assembly=b37>
##contig=<ID=GL000213.1,length=164239,assembly=b37>
##contig=<ID=GL000211.1,length=166566,assembly=b37>
##contig=<ID=GL000199.1,length=169874,assembly=b37>
##contig=<ID=GL000217.1,length=172149,assembly=b37>
##contig=<ID=GL000216.1,length=172294,assembly=b37>
##contig=<ID=GL000215.1,length=172545,assembly=b37>
##contig=<ID=GL000205.1,length=174588,assembly=b37>
##contig=<ID=GL000219.1,length=179198,assembly=b37>
##contig=<ID=GL000224.1,length=179693,assembly=b37>
##contig=<ID=GL000223.1,length=180455,assembly=b37>
##contig=<ID=GL000195.1,length=182896,assembly=b37>
##contig=<ID=GL000212.1,length=186858,assembly=b37>
##contig=<ID=GL000222.1,length=186861,assembly=b37>
##contig=<ID=GL000200.1,length=187035,assembly=b37>
##contig=<ID=GL000193.1,length=189789,assembly=b37>
##contig=<ID=GL000194.1,length=191469,assembly=b37>
##contig=<ID=GL000225.1,length=211173,assembly=b37>
##contig=<ID=GL000192.1,length=547496,assembly=b37>
##contig=<ID=NC_007605,length=171823,assembly=b37>
##contig=<ID=hs37d5,length=35477943,assembly=b37>
