#!/usr/bin/perl
#============================================================================
# Name        		: surroundingSequence.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Thu Sep 29 14:51:22 2011
# Last Modified By	: 
# Last Modified On	: Thu Sep 29 14:51:22 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl surroundingSequence.pl [option] <fasta> <infile>

	-t	vcf or bed [default vcf]
	-s	surrounding size [defualt 10]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Bio::DB::Sam;
use Data::Dumper;

my ($in,$out);
my ($opt_h,$opt_t,$opt_s);
GetOptions("h"	=>\$opt_h,"t=s"=>\$opt_t,"s=i"=>\$opt_s);
if(@ARGV<2 || $opt_h) { usage(); }
my $fafile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_t)) { $opt_t="vcf"; }
if($opt_t ne "vcf" && $opt_t ne "bed") { usage(); }
if(!defined($opt_s)) { $opt_s=10; }

if($infile=~/\.gz$/) { open $in,"bgzip -cd | $infile" or die "Cann't open file $infile ($!)\n"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
my $fai = Bio::DB::Sam::Fai->load($fafile);

my %context=('all'=>[],"A"=>[],"T"=>[],"C"=>[],"G"=>[]);
for(my $i=0;$i<$opt_s*2+1;$i++)
{
	push @{$context{'all'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'A'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'T'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'C'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'G'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
}
open $in,$infile or die "Cann't open file $infile ($!) \n";
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
	my $dna_string = $fai->fetch(sprintf("$chr:%d-%d",$pos-$opt_s,$pos+$opt_s));
	$dna_string="\U$dna_string";
	my @ss=split //,$dna_string;
	for(my $i=0;$i<@ss;$i++)
	{
		$context{'all'}->[$i]->{$ss[$i]}++;
		$context{$ref}->[$i]->{$ss[$i]}++;
	}
	#printf "$chr\t%d\t%d\t%s\n",$pos-1,$pos,$dna_string;
}
foreach my $k ("A","T","C","G","all")
{
	printf "mutation context of $k:";
	for(my $_ii=-$opt_s;$_ii<=$opt_s;$_ii++) { printf "\t%d",$_ii; }
	printf "\n";
	my $p=$context{$k};
	foreach my $i ("A","T","C","G")
	{
		printf "$i";
		for(my $j=0;$j<@$p;$j++)
		{
			printf "\t%s",$p->[$j]->{$i};
		}
		printf "\n";
	}
}

############################################################################
sub usage
{
	die `pod2text $0`;
}
