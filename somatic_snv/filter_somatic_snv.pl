#!/usr/bin/env perl
#============================================================================
# Name        		: filter_somatic_snv.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Thu Feb 28 23:07:58 2013
# Last Modified By	: 
# Last Modified On	: Thu Feb 28 23:07:58 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl filter_somatic_snv.pl [option] 

	-CLR		CLR less than this will be discarded [default 15]
	-normalDP	normal depth less than this will be discarded [default 10]
	-tumorDP	tumor depth less than this will be discarded [default 10]
	-SB		strand bias score larger than this will be discarded [default 40]
	-v		verbose [default: OFF]
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_CLR,$opt_normalDP,$opt_tumorDP,$opt_SB,$opt_v);
GetOptions("h"	=>\$opt_h,"CLR=f"=>\$opt_CLR,"normalDP=f"=>\$opt_normalDP,"tumorDP=f"=>\$opt_tumorDP,"SB=f"=>\$opt_SB,"v"=>\$opt_v);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_CLR)) { $opt_CLR=15; }
if(!defined($opt_normalDP)) { $opt_normalDP=10; }
if(!defined($opt_tumorDP)) { $opt_tumorDP=10; }
if(!defined($opt_SB)) { $opt_SB=40; }

while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { print "$_\n"; next; }
	#chr1    11851003        .       G       C       153     .       DP=312;VDB=0.0120;AF1=0.25;CLR=186;AC1=1;DP4=152,70,55,2;MQ=47;FQ=154;PV4=2.3e-06,0.33,1.4e-13,1        GT:PL:DP:SP:GQ  0/0:0,255,255:167:0:99  0/1:186,0,255:112:40:99
	my @F=split /\t/,$_;
	my @tag=split /:/,$F[-3];
	my %tag=();
	for(my $i=0;$i<@tag;$i++)
	{
		$tag{$tag[$i]}=$i;
	}
	my @N=split /:/,$F[-2];
	my @T=split /:/,$F[-1];
	my $normalGT=$N[$tag{'GT'}];
	my $tumorGT=$T[$tag{'GT'}];
	my $normalDP=$N[$tag{'DP'}];
	my $tumorDP=$T[$tag{'DP'}];
	my $normalSB=$N[$tag{'SP'}];
	my $tumorSB=$T[$tag{'SP'}];
	my ($CLR) = /CLR=(.+?);/;
	
	if(!defined($CLR)) 
	{
		if($opt_v) { printf STDERR "##missing data in one sample:\t$line\n"; }
		next; 
	}
	if(!($normalGT eq "0/0" && ($tumorGT eq "0/1" || $tumorGT eq "1/1")))
	{ 
		if($opt_v) { printf STDERR "##GT:\t$line\n"; }
		next; 
	}
	if($CLR<$opt_CLR) 
	{ 
		if($opt_v) { printf STDERR "##CLR:\t$line\n"; }
		next; 
	}
	if($normalDP<$opt_normalDP) 
	{ 
		if($opt_v) { printf STDERR "##normal depth:\t$line\n"; }
		next; 
	}
	if($tumorDP<$opt_tumorDP) 
	{ 
		if($opt_v) { printf STDERR "##tumor depth:\t$line\n"; }
		next; 
	}
	if($normalSB<$opt_SB && $tumorSB>$opt_SB) 
	{ 
		if($opt_v) { printf STDERR "##strand bias:\t$line\n"; }
		next; 
	}
	
	printf "$line\n";

}



############################################################################
sub usage
{
	die `pod2text $0`;
}
