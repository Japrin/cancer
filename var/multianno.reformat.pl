#!/usr/bin/perl
#============================================================================
# Name        		: multianno.reformat.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Mon Dec 23 16:04:25 2013
# Last Modified By	: 
# Last Modified On	: Mon Dec 23 16:04:25 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl multianno.reformat.pl [option] <infile>

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

my $title_line="";
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
	#CopyNumber=3;Size=1735587;SVID=10;SVType=gain
	#Only the first two lines
	chomp;
	$title_line=$_;
	$_=<$in>;
	chomp;
	my @F=split /\t/;
	my @aa=split /[=;]/,$F[-1];
	my @bb=();
	for(my $i=0;$i<@aa;$i=$i+2)
	{
		push @bb,$aa[$i];
	}
	my $_b=join("\t",@bb);
	$title_line=~s/Otherinfo/$_b/;
	last;
}
close $in;
####
print "$title_line\n";
open $in,$infile or die "Cann't open file $infile ($!) \n";
<$in>;
while(<$in>)
{
	chomp;
	my @F=split /\t/;
	my @aa=split /[=;]/,$F[-1];
	my @bb=();
	for(my $i=1;$i<@aa;$i=$i+2)
	{
		push @bb,$aa[$i];
	}
	my $_b=join("\t",@bb);
	$F[-1]=$_b;
	print join("\t",@F)."\n";
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
