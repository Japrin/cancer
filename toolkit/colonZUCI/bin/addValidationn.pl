#!/usr/bin/perl
#============================================================================
# Name        		: addValidationn.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Thu Sep  8 12:30:28 2011
# Last Modified By	: 
# Last Modified On	: Thu Sep  8 12:30:28 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl addValidationn.pl [option] <vcffile> <valfile>

	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<2 || $opt_h) { usage(); }
my $infile=shift @ARGV;
my $valfile=shift @ARGV;

my %list=();
readList(\%list,$valfile);

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
			print "##INFO=<ID=validation,Number=.,Type=String,Description=\"validation status(Yes,true somatic; No,false somatic; NA, not validated\">\n";
			print "$_\n";
			next;
		}
		if(/^##INFO/) { push @infoHeader,$_; next; }
		else { print "$_\n"; next; }
	}
	##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Tumor
	my @field=split /\t/;
	my ($chr,$pos)=@field[0,1];
	$field[7]=~s/;$//;
	if(exists($list{"$chr:$pos"})) { $field[7].=sprintf(";validation=%s;",$list{"$chr:$pos"}->{'somatic'}); }
	else { $field[7].=sprintf(";validation=NA;"); }
	printf "%s\n",join("\t",@field);
}



############################################################################
sub readList
{
	my ($list,$infile)=@_;
	my ($in);
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		##Chr    Start   Ref     Obs     somatic
		my @field=split /\t/;
		my ($chr,$pos,$ref,$alt,$somatic)=@field[0,1,2,3,4];
		my $key="$chr:$pos";
		$list->{$key}={'ref'=>$ref,'alt'=>$alt,'somatic'=>$somatic};
	}
}
sub usage
{
	die `pod2text $0`;
}
