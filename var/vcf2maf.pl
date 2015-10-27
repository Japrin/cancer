#!/usr/bin/env perl
#============================================================================
# Name        		: vcf2maf.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Tue May 13 19:38:15 2014
# Last Modified By	: 
# Last Modified On	: Tue May 13 19:38:15 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl vcf2maf.pl [option] <infile> [outfile]

	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;
#my $outfile=shift @ARGV;

if($infile=~/.gz$/) { open $in,"gzip -cd $infile |" or die "$!"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
