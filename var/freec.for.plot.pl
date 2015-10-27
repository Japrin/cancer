#!/usr/bin/env perl
#============================================================================
# Name        		: freec.for.plot.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Mon Dec 23 17:23:54 2013
# Last Modified By	: 
# Last Modified On	: Mon Dec 23 17:23:54 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl freec.for.plot.pl [option] <infile>

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
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;

###determin bin width
open $in,$infile or die "Cann't open file $infile ($!) \n";
$_=<$in>;	## header line
$_=<$in>;	
chomp;
my $line1=$_;	## line 1
$_=<$in>;
chomp;
my $line2=$_;	## line 2
close $in;
my $beg1=(split /\t/,$line1)[1];
my $beg2=(split /\t/,$line2)[1];
my $bin_width=$beg2-$beg1;
###
open $in,$infile or die "Cann't open file $infile ($!) \n";
$_=<$in>;
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#Chromosome      Start   Ratio   MedianRatio     CopyNumber
	#1       1       -1      -1      2
	#1       64282   1.77735 1.10418 2
	my ($chr,$beg,$bin_ratio,$seg_ratio,$cn)=@F;
	printf "$chr\t$beg\t%d\t%s\n",$beg+$bin_width-1,$bin_ratio==-1?"NA":log($bin_ratio)/log(2);
	printf STDERR "$chr\t$beg\t%d\t%s\n",$beg+$bin_width-1,$seg_ratio==-1?"NA":log($seg_ratio)/log(2);

}
###################################################################




sub usage
{
	die `pod2text $0`;
}
