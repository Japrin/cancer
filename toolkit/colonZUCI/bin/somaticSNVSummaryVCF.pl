#!/usr/bin/perl
#============================================================================
# Name        		: somaticSNVSummaryVCF.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Wed Sep  7 17:01:47 2011
# Last Modified By	: 
# Last Modified On	: Wed Sep  7 17:01:47 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl somaticSNVSummaryVCF.pl [option] <infile>

	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;

my %precedenceFunc=('exonic'=>1,
		'exonic,splicing'=>1,
		'splicing'=>2,
		'ncRNA'=>3,
		'ncRNA_exonic'=>3.1,
		'ncRNA_intronic'=>3.2,
		'ncRNA_UTR5'=>3.3,
		'ncRNA_UTR3'=>3.4,
		'UTR5'=>4,
		'UTR3'=>5,
		'intronic'=>6,
		'upstream'=>7,
		'downstream'=>8,
		'intergenic'=>9,
		);
my %precedenceExonicFunc=(
		'stopgain SNV'=>1,
		'stoploss SNV'=>2,
		'nonsynonymous SNV'=>3,
		'synonymous SNV'=>4,
		);

my %list=();
if($infile =~ /\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!)\n"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
while(<$in>)
{
	chomp;
	if(/^\s*$/ || /^#/) { next; }
	if(/\bExonicFunc=(.+?);/) { $list{'ExonicFunc'}->{$1}++; }
	if(/\bFunc=(.+?);/) { $list{'Func'}->{$1}++; }
}
print "------Func------\n";
foreach (sort { $precedenceFunc{$a}<=>$precedenceFunc{$b} } keys %{$list{'Func'}})
{
	printf "$_\t%d\n",$list{'Func'}->{$_};
}
print "------ExonicFunc------\n";
foreach (sort { $precedenceExonicFunc{$a}<=>$precedenceExonicFunc{$b} } keys %{$list{'ExonicFunc'}})
{
	printf "$_\t%d\n",$list{'ExonicFunc'}->{$_};
}

############################################################################
sub usage
{
	die `pod2text $0`;
}
