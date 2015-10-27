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

	perl SNVSpectrumVCF.pl [option] <vcffile>

	-f=<dir>	reference data dir
	-c=<int>	output spectrum mode:14, 7, 9.[default 7]
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use lib "/home/zhoudonger/01bin/lib/perl/lib/perl/5.10.1";
use Bio::DB::Sam;
use Data::Dumper;

my ($in,$out);
my ($opt_h,$opt_f,$opt_c);
GetOptions("h"	=>\$opt_h,"f=s"=>\$opt_f,"c=i"=>\$opt_c);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;
if($opt_f && ! -e $opt_f) { usage(); }
if(!defined($opt_c)) { $opt_c=7; }
if($opt_c !=14 && $opt_c !=7 && $opt_c !=9) { usage(); }

my %spectrum=(
	"A->T"=>0,
	"T->A"=>0,

	"A->C"=>0,
	"T->G"=>0,

	"A->G"=>0,
	"T->C"=>0,

	"C->A"=>0,
	"G->T"=>0,
	"C->AinCpG"=>0,
	"G->TinCpG"=>0,
	"C->AnotCpG"=>0,
	"G->TnotCpG"=>0,

	"C->T"=>0,
	"G->A"=>0,
	"C->TinCpG"=>0,
	"C->TnotCpG"=>0,
	"G->AinCpG"=>0,
	"G->AnotCpG"=>0,

	"C->G"=>0,
	"G->C"=>0,
	"C->GinCpG"=>0,
	"G->CinCpG"=>0,
	"C->GnotCpG"=>0,
	"G->CnotCpG"=>0,
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
	my ($chr,$pos,$ref,$alt)=@field[0,1,3,4];
	my $dna_string = $fai->fetch(sprintf("$chr:%d-%d",$pos-1,$pos+1));
	$dna_string="\U$dna_string";
	my ($preBase,$curBase,$nextBase)=split //,$dna_string;
	if($curBase ne $ref) { printf STDERR "inconsistent ref\t$_\n"; next; }

	my $key=sprintf("%s->%s",$curBase,$alt);
	$spectrum{$key}++;
	if($curBase eq "C" && $nextBase eq "G")
	{
		$key=sprintf("%s->%sinCpG",$curBase,$alt);
		$spectrum{$key}++;
	}elsif($preBase eq "C" && $curBase eq "G")
	{
		$key=sprintf("%s->%sinCpG",$curBase,$alt);
		$spectrum{$key}++;
	}elsif($curBase eq "C" or $curBase eq "G")
	{
		$key=sprintf("%s->%snotCpG",$curBase,$alt);
		$spectrum{$key}++;
	}
}
if($opt_c == 7)
{
	printf "A:T->T:A\t%d\n",$spectrum{"A->T"}+$spectrum{"T->A"};
	printf "A:T->C:G\t%d\n",$spectrum{"A->C"}+$spectrum{"T->G"};
	printf "A:T->G:C\t%d\n",$spectrum{"A->G"}+$spectrum{"T->C"};
	printf "C:G->A:T\t%d\n",$spectrum{"C->A"}+$spectrum{"G->T"};
	printf "C:G in CpG -> T:A\t%d\n",$spectrum{"C->TinCpG"}+$spectrum{"G->AinCpG"};
	printf "C:G not CpG -> T:A\t%d\n",$spectrum{"C->TnotCpG"}+$spectrum{"G->AnotCpG"};
	printf "C:G->G:C\t%d\n",$spectrum{"C->G"}+$spectrum{"G->C"};
}elsif($opt_c == 14)
{
	printf "A->T\t%d\n",$spectrum{"A->T"};
	printf "T->A\t%d\n",$spectrum{"T->A"};
	printf "A->C\t%d\n",$spectrum{"A->C"};
	printf "T->G\t%d\n",$spectrum{"T->G"};
	printf "A->G\t%d\n",$spectrum{"A->G"};
	printf "T->C\t%d\n",$spectrum{"T->C"};
	printf "C->A\t%d\n",$spectrum{"C->A"};
	printf "G->T\t%d\n",$spectrum{"G->T"};
	printf "C->TinCpG\t%d\n",$spectrum{"C->TinCpG"};
	printf "G->AinCpG\t%d\n",$spectrum{"G->AinCpG"};
	printf "C->TnotCpG\t%d\n",$spectrum{"C->TnotCpG"};
	printf "G->AnotCpG\t%d\n",$spectrum{"G->AnotCpG"};
	printf "C->G\t%d\n",$spectrum{"C->G"};
	printf "G->C\t%d\n",$spectrum{"G->C"};
}elsif($opt_c == 9)
{
	printf "A:T->T:A\t%d\n",$spectrum{"A->T"}+$spectrum{"T->A"};
	printf "A:T->C:G\t%d\n",$spectrum{"A->C"}+$spectrum{"T->G"};
	printf "A:T->G:C\t%d\n",$spectrum{"A->G"}+$spectrum{"T->C"};
	printf "C:G in CpG -> A:T\t%d\n",$spectrum{"C->AinCpG"}+$spectrum{"G->TinCpG"};
	printf "C:G not CpG -> A:T\t%d\n",$spectrum{"C->AnotCpG"}+$spectrum{"G->TnotCpG"};
	printf "C:G in CpG -> T:A\t%d\n",$spectrum{"C->TinCpG"}+$spectrum{"G->AinCpG"};
	printf "C:G not CpG -> T:A\t%d\n",$spectrum{"C->TnotCpG"}+$spectrum{"G->AnotCpG"};
	printf "C:G in CpG -> G:C\t%d\n",$spectrum{"C->GinCpG"}+$spectrum{"G->CinCpG"};
	printf "C:G not CpG -> G:C\t%d\n",$spectrum{"C->GnotCpG"}+$spectrum{"G->CnotCpG"};
}
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
