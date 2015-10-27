#!/usr/bin/env perl
#============================================================================
# Name        		: freec.loh.region.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Fri Jan 10 00:26:32 2014
# Last Modified By	: 
# Last Modified On	: Fri Jan 10 00:26:32 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl freec.loh.region.pl [option] <infile>
	
	-f	loh threshold [default 0.6]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_f);
GetOptions("h"	=>\$opt_h,"f=f"=>\$opt_f);
if(@ARGV<0 || $opt_h) { usage(); }
#my $infile=shift @ARGV;
if(!defined($opt_f)) { $opt_f=0.6; }

my $preChr="";
my $preBeg="";
my $preEnd="";
my $prePos="";
my $preA="";
my $preB="";
my $preInfo="";
#open $in,$infile or die "Cann't open file $infile ($!) \n";
my $h=<>;
#print "$h";
while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#Chromosome      Position        BAF     FittedA FittedB A       B       uncertainty
	#1       755384  0.0645161       0.11    0.89    0       1       0.0363972
	#1       16931554        0.654412        0.35047 0.64953 0.333333        0.666667        4.80007 1       16931376        17235504
	#        6       gain    AAAABB  NA      somatic NA
	my ($chr,$pos,$A,$B)=@F[0,1,5,6];
	my $info=join("\t",@F[8..16]);
	if($A eq "-1") 
	{
		if($preChr ne "" && $preA ne "-1")
		{
			##output
			#if($preA>=$opt_f || $preB>=$opt_f)
			if($preA != 0.5)
			{
				print "$preChr\t$preBeg\t$preEnd\t$preA\t$preB\t$preInfo\n";
			}
		}
		$preChr=$chr;
		$preBeg=$pos;
		$preEnd=$pos;
		$prePos=$pos;
		$preA=$A;
		$preB=$B;
		$preInfo=$info;
		next;
	}

	if($chr eq $preChr && $A eq $preA && $B eq $preB)
	{
		$prePos=$pos;
		$preEnd=$prePos;
	}else
	{
		if($preChr ne "")
		{
			##output
			#print "$preChr\t$preBeg\t$preEnd\t$preA\t$preB\n";
		}
		$preChr=$chr;
		$preBeg=$pos;
		$preEnd=$pos;
		$prePos=$pos;
		$preA=$A;
		$preB=$B;
		$preInfo=$info;
	}
}
if($preChr ne "" && $preA ne "-1")
{
	##output
	#if($preA>=$opt_f || $preB>=$opt_f)
	if($preA != 0.5)
	{
		print "$preChr\t$preBeg\t$preEnd\t$preA\t$preB\t$preInfo\n";
	}
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
