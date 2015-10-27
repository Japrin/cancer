#!/usr/bin/perl
#============================================================================
# Name        		: filter_snp_GATK.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Sat Aug  6 17:36:54 2011
# Last Modified By	: 
# Last Modified On	: Sat Aug  6 17:36:54 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl filter_snp_GATK.pl [option] <vcffile>

	-Q		Alt quailty threshold [default 20]
	-MQ		RMS Mapping quality threshold [default 10]
	-MQ0		Number of MQ0 reads threshold [default 8]
	--snv		filter out those in dbSNP or 1KG
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out,$out_l);
my ($opt_h,$QTh,$MQTh,$MQ0Th,$opt_N,$opt_l,$opt_snv);
GetOptions("h"	=>\$opt_h,"Q=i"=>\$QTh,"MQ=i"=>\$MQTh,"MQ0=i"=>\$MQ0Th,"snv"=>\$opt_snv);
if(@ARGV<1 || $opt_h) { usage(); }
my $vcffile=shift @ARGV;
if(!defined($QTh)) { $QTh=20; }
if(!defined($MQTh)) { $MQTh=10; }
if(!defined($MQ0Th)) { $MQ0Th=8; }

if($vcffile =~ /\.gz$/) { open $in,"bgzip -cd $vcffile|" or die "Cann't open file $vcffile ($!) \n"; }
else { open $in,$vcffile or die "Cann't open file $vcffile ($!) \n"; }
while(<$in>)
{
	chomp;
	if(/^\s*$/ || /^#/) 
	{
		print "$_\n";
		next; 
	}
	##chr1    120     .       T       C       12.35   HARD_TO_VALIDATE;LowQual        Func=intergenic;Gene=NONE(dist=NONE),LOC100288778(dist=4105);AC=1;AF=0.50;AN=2;BaseQRankSum=0.423;DP=20;Dels=0.00;HRun=3;HaplotypeScore=2.9841;MQ=18.06;MQ0=6;MQRankSum=0.085;QD=0.62;ReadPosRankSum=1.944;SB=2.38      GT:AD:DP:GQ:PL  0/1:16,3:20:42.08:42,0,182
	my @field=split /\t/;
	my ($chr,$pos,$rsID,$ref,$alt,$qual,$filter)=@field[0,1,2,3,4,5,6];
	if($rsID ne "." && $opt_snv) { next; }
	if(/1000G.*?=(.+?);/ && $opt_snv) { next; }
	my ($_MQ)= /MQ=(.+?);/;
	if($qual<$QTh) { next; }
	if($_MQ<$MQTh) { next; }
	if(/MQ0=(.+?);/) { if($1>=$MQ0Th){ next; } }
	if(length($ref)!=length($alt)) { next; }
	
	#if($filter eq "." || $filter eq "PASS" || $filter eq "TruthSensitivityTranche99.00to99.90")
	if($filter eq "." || $filter eq "PASS")
	{
		print "$_\n";
	}else
	{
		#printf STDERR "Filtered SNP:\t$_\n";
	}
}

sub usage
{
	die `pod2text $0`;
}
