#!/usr/bin/perl
#============================================================================
# Name        		: addContext.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Thu Sep 29 18:47:03 2011
# Last Modified By	: 
# Last Modified On	: Thu Sep 29 18:47:03 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl addContext.pl [option] <fastaFile> <infile>

	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Bio::DB::Sam;

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<2 || $opt_h) { usage(); }
my $fafile=shift @ARGV;
my $infile=shift @ARGV;

if($infile=~/\.gz$/) { open $in,"bgzip -cd | $infile" or die "Cann't open file $infile ($!)\n"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
my $fai = Bio::DB::Sam::Fai->load($fafile);

my @infoHeader=();
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
	chomp;
	if(/^#/) 
	{ 
		if(/^#CHROM/)
		{
			printf "%s\n",join("\n",@infoHeader);
			print "##INFO=<ID=context1,Number=.,Type=String,Description=\"context1(NpC)\">\n";
			print "##INFO=<ID=context2,Number=.,Type=String,Description=\"context2(CpN)\">\n";
			print "$_\n";
			next;
		}
		if(/^##INFO/) { push @infoHeader,$_; next; }
		else { print "$_\n"; next; }
	}
	##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Tumor
	my @field=split /\t/;
	my ($chr,$pos,$ref,$alt);
	($chr,$pos,$ref,$alt)=@field[0,1,3,4];
	$field[7]=~s/;$//;
	my $dna_string = $fai->fetch(sprintf("$chr:%d-%d",$pos-1,$pos+1));
	$dna_string="\U$dna_string";
	my @ss=split //,$dna_string;
	if($ss[1] ne $ref)
	{
		warn "ref($ref) not match fast($ss[1])\t$_\n";
		next;
	}
	$field[7].=sprintf(";context1=%s;context2=%s","$ss[0]$ss[1]","$ss[1]$ss[2]");
	printf "%s\n",join("\t",@field);
}


############################################################################
sub usage
{
	die `pod2text $0`;
}
