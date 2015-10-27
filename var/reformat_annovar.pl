#!/usr/bin/perl -w 

=pod

=head1 NAME

			reformat_annovar.pl

=head1 SYNOPSIS

perl reformat_annovar.pl <infile> -id [ID]

=head1 DESCRIPTION

reformat table_annovar output file to pretty table format

=head1 OPTIONS

=over

=item -id

specify the ID of the samples, must be in the same order with annovar input

=back

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($id,$h,$m); 
GetOptions(
	"id=s"=>\$id,
	"h|help"=>\$h,
	"m|man"=>\$m
);

if ($h) {pod2usage(-verbose => 1);}
if ($m) {pod2usage(-verbose  => 2);}
if (!defined $h and !defined $m and !defined $id)  { pod2usage(-verbose=>0);}


my $infile = shift;
my @id = split ",", $id;

open In, "$infile" or die "$!";
my @fields = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER");
my ($length,%header); 

while (<In>) {
	chomp;
	if (/^Chr/) {
		s/\.knownGene//g;
		my @new_fields = split "\t", $_;
		$length = @new_fields;
		splice (@new_fields, 0, 5);
		pop @new_fields;
		@fields = (@fields, @new_fields,"INFO","FORMAT",@id);
		print join("\t",@fields)."\n";
	}
	else {
		my @records = split "\t", $_;
		## "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"
		print join("\t",@records[($length)..($length+6)]);
		# not necessary now:
		#if ($records[5] eq "splicing") {
		#	if ($records[6] =~ /\((.+?)\)/) {
		#		$records[8] = $1;
		#		$records[6] =~ s/\(.+?\)//;
		#	}
		#}
		## Func.knownGene  Gene.knownGene  GeneDetail.knownGene    ExonicFunc.knownGene    AAChange.knownGene      cpgIslandExt    evofold wgRna   targetScanS     phastConsElements46way  genomicSuperDups        tfbsConsSites   cytoband	esp6500si_all   1000g2012apr_all        snp138  ljb23_sift      ljb23_pp2hdiv   ljb23_pp2hvar   gerp++elem      cosmic68
		print "\t".join("\t",@records[5..($length-2)]);	
		#foreach (5..$length-2) {
		#	print "\t$records[$_]";
		#}
		# INFO, FORMAT and SAMPLE:
		print "\t$records[$length+7]";
		foreach (($length+8)..$#records) {
			print "\t$records[$_]";
		}
		print "\n";
	}
}
close In;
