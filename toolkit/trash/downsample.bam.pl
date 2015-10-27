#!/usr/bin/perl
#============================================================================
# Name        		: downsample.bam.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Thu Apr 25 15:16:08 2013
# Last Modified By	: 
# Last Modified On	: Thu Apr 25 15:16:08 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl downsample.bam.pl [option] [infile]

	-f	only chromosomes: chr1,...chr22,chrX, chrY [default all]
	-a	downsample to 1/a [default 30]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_f,$opt_a);
GetOptions("h"	=>\$opt_h,"f"=>\$opt_f,"a=i"=>\$opt_a);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_a)) { $opt_a=30; }

my %CHR_LIST=("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y");

while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	my ($flag,$chr)=@F[1,2];
	$chr=~s/^chr//;
	if($opt_f)
	{
		if(!exists($CHR_LIST{$chr})) { next; }
	}
	# HWI-ST904:164:C0K5RACXX:4:1305:15328:24030      pPR1    chrM    1       60      101M    =       171     271     GATC
	my $r=int(rand($opt_a));
	if($r==0)
	{
		print "$line\n";
	}
}
###################################################################




sub usage
{
	die `pod2text $0`;
}

