#!/usr/bin/perl
#============================================================================
# Name        		: SNVSpectrumVCF.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Fri Sep  9 12:32:47 2011
# Last Modified By	: 
# Last Modified On	: Fri Sep  9 12:32:47 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl SNVSpectrumVCF.pl [option] <snvfile>

	-f=<file>	reference data dir
	-t=<string>	file type, vcf or bed [default vcf]
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use lib "/ifshk1/BC_CANCER/01bin/lib/perl/local/lib/perl5/site_perl";
use lib "/ifshk1/BC_CANCER/01bin/lib/perl/local/lib/perl5/x86_64-linux-thread-multi/"; 
use Bio::DB::Sam;
use Data::Dumper;

my ($in,$out);
my ($opt_h,$opt_f,$opt_c,$opt_t);
GetOptions("h"	=>\$opt_h,"f=s"=>\$opt_f,"c=i"=>\$opt_c,"t=s"=>\$opt_t);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;
if($opt_f && ! -e $opt_f) { usage(); }
if(!defined($opt_c)) { $opt_c=7; }
if(!defined($opt_t)) { $opt_t="vcf"; }
if($opt_c !=14 && $opt_c !=7) { usage(); }
if($opt_t ne "vcf" && $opt_t ne "bed") { usage(); }

my %spectrum=(
	"A->T"=>0,
	"T->A"=>0,
	"A->C"=>0,
	"T->G"=>0,
	"A->G"=>0,
	"T->C"=>0,
	"C->A"=>0,
	"G->T"=>0,
	"C->G"=>0,
	"G->C"=>0,
	"CpG"=>0,
	"TpC"=>0,
	"GpA"=>0,
	"Total"=>0,
);
	#A:T->T:A	A->T || T->A
	#A:T->C:G	A->C || T->G
	#A:T->G:C	A->G || T->C
	#C:G->A:T	C->A || G->T
	#C:G->T:A	C->T || G->A	(C->T in CpG || G->A in CpG) (not in CpG)
	#C:G->G:C	C->G || G->C

if($infile=~/\.gz$/) { open $in,"bgzip -cd | $infile" or die "Cann't open file $infile ($!)\n"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
my $fai = Bio::DB::Sam::Fai->load($opt_f);
while(<$in>)
{
	chomp;
	if(/^\s*$/ || /^#/) { next; }
	##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Adenoma
	my @field=split /\t/;
	my ($chr,$pos,$ref,$alt);
	if($opt_t eq "vcf")
	{
		($chr,$pos,$ref,$alt)=@field[0,1,3,4];
	}elsif($opt_t eq "bed")
	{
		($chr,$pos,$ref,$alt)=@field[0,2,3,4];
	}
	my $dna_string = $fai->fetch(sprintf("$chr:%d-%d",$pos-1,$pos+1));
	$dna_string="\U$dna_string";
	my ($preBase,$curBase,$nextBase)=split //,$dna_string;
	if($curBase ne $ref) { printf STDERR "inconsistent ref\t$_\n"; next; }
	my $key=sprintf("%s->%s",$curBase,$alt);
	$spectrum{$key}++;
	if($curBase eq "C" && $nextBase eq "G")
	{
		$spectrum{"CpG"}++;
	}elsif($preBase eq "C" && $curBase eq "G")
	{
		$spectrum{"CpG"}++;
	}elsif($preBase eq "T" && $curBase eq "C")
	{
		$spectrum{"TpC"}++;
	}elsif($curBase eq "G" && $nextBase eq "A")
	{
		$spectrum{"GpA"}++;
	}
	$spectrum{"Total"}++;
}
printf "A:T->T:A\t%d\n",$spectrum{"A->T"}+$spectrum{"T->A"};
printf "A:T->C:G\t%d\n",$spectrum{"A->C"}+$spectrum{"T->G"};
printf "A:T->G:C\t%d\n",$spectrum{"A->G"}+$spectrum{"T->C"};
printf "C:G->A:T\t%d\n",$spectrum{"C->A"}+$spectrum{"G->T"};
printf "C:G->T:A\t%d\n",$spectrum{"C->T"}+$spectrum{"G->A"};
printf "C:G->G:C\t%d\n",$spectrum{"C->G"}+$spectrum{"G->C"};
printf "CpG\t%d\n",$spectrum{"CpG"};
printf "TpC\t%d\n",$spectrum{"TpC"}+$spectrum{"GpA"};
printf "Total\t%d\n",$spectrum{"Total"};
	
	#A:T->T:A	A->T || T->A
	#A:T->C:G	A->C || T->G
	#A:T->G:C	A->G || T->C
	#C:G->A:T	C->A || G->T
	#C:G->T:A	C->T || G->A	(C->T in CpG || G->A in CpG) (not in CpG)
	#C:G->G:C	C->G || G->C

#print Dumper(%spectrum);

############################################################################
sub usage
{
	die `pod2text $0`;
}
