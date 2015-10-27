#!/usr/bin/perl
#============================================================================
# Name        		: MentoCaloCalRatio.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Wed Sep 14 22:48:04 2011
# Last Modified By	: 
# Last Modified On	: Wed Sep 14 22:48:04 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl MentoCaloCalRatio.pl [option] <variant file> <exonic variant file>

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
my $variantFile=shift @ARGV;
my $exonicVariantFile=shift @ARGV;

my %res=();
open $in,$variantFile or die "Cann't open file $variantFile ($!) \n";
while(<$in>)
{
	chomp;
	if(/^\s*$/) { next; }
	#exonic  SAMD11  chr1    856305  856305  A       G       4726    A.bed.gz        781
	my @field=split /\t/;
	my ($func,$iSim)=@field[0,7];
	exists($res{$iSim}) or $res{$iSim}={"NonSyn"=>0,"Syn"=>0};
	if($func eq "splicing") { $res{$iSim}->{"NonSyn"}++; }
}
open $in,$exonicVariantFile or die "Cann't open file $exonicVariantFile ($!) \n";
while(<$in>)
{
	chomp;
	if(/^\s*$/) { next; }
	#line5   nonsynonymous SNV       SAMD11:NM_152486:exon12:c.A1624C:p.K542Q,       chr1    868555  868555  A       C       1964    A.bed.gz        999
	my @field=split /\t/;
	my ($exonic_func,$iSim)=@field[1,8];
	exists($res{$iSim}) or $res{$iSim}={"NonSyn"=>0,"Syn"=>0};
	if($exonic_func eq "synonymous SNV") { $res{$iSim}->{"Syn"}++; }
	else{ $res{$iSim}->{"NonSyn"}++; }
}
foreach (keys %res)
{
	printf "$_\t%4.4f\n",$res{$_}->{"Syn"}?$res{$_}->{"NonSyn"}/$res{$_}->{"Syn"}:"NA";
}


############################################################################
sub usage
{
	die `pod2text $0`;
}
