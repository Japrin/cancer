#!/usr/bin/perl
#============================================================================
# Name        		: GATKSomaticIndel_filter.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Wed Aug 10 17:12:50 2011
# Last Modified By	: 
# Last Modified On	: Wed Aug 10 17:12:50 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl GATKSomaticIndel_filter.pl [option] <input vcf>

	--RStartTh	Threshold for RStart [default 5]
	--REndTh	Threshold for REnd [default 5]
	--MMTh		Threshold for MM [default 1]
	--MQSTh		Threshold for MQS [default 30]
	--NQSBQTh	Threshold for NQSBQ [default 20]
	--NQSMMTh	Threshold for NQSMM [default 0.01]
	--NDepthTh	minimum depth required for normal [default 10]
	--TDepthTh	minimum depth required for tumor [default 10]
	--NFreqTh	maximum variant frequency for normal [default 0.05]
	--TFreqTh	minimum variant frequency for tumor [default 0.35]
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$RStartTh,$REndTh,$MMTh,$MQSTh,$NQSBQTh,$NQSMMTh,$NDepthTh,$TDepthTh,$NFreqTh,$TFreqTh,$opt_NIndex,$opt_TIndex);
GetOptions("h"	=>\$opt_h,"RStartTh=i"=>\$RStartTh,"REndTh=i"=>\$REndTh,"MMTh=f"=>\$MMTh,"MQSTh=f"=>\$MQSTh,"NQSBQTh=f"=>\$NQSBQTh,"NQSMMTh=f"=>\$NQSMMTh,"NDepthTh=i"=>\$NDepthTh,"TDepthTh=i"=>\$TDepthTh,"NFreqTh=f"=>\$NFreqTh,"TFreqTh=f"=>\$TFreqTh,"NIndex=i"=>\$opt_NIndex,"TIndex=i"=>\$opt_TIndex);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;
if(!defined($RStartTh)) { $RStartTh=5; }
if(!defined($REndTh)) { $REndTh=5; }
if(!defined($MMTh)) { $MMTh=1; }
if(!defined($MQSTh)) { $MQSTh=30; }
if(!defined($NQSBQTh)) { $NQSBQTh=20; }
if(!defined($NQSMMTh)) { $NQSMMTh=0.01; }
if(!defined($NDepthTh)) { $NDepthTh=10; }
if(!defined($TDepthTh)) { $TDepthTh=10; }
if(!defined($NFreqTh)) { $NFreqTh=0.05; }
if(!defined($TFreqTh)) { $TFreqTh=0.35; }
#if(!defined($opt_NIndex)) { $opt_NIndex=-2; }
#if(!defined($opt_TIndex)) { $opt_TIndex=-1; }

my $auto=0;
my @sampleName=();

if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "$!"; }
else { open $in,$infile or die "$!"; }
while(<$in>)
{
	chomp;
	my $line=$_;
	my @F=split /\t/;
	#GT:AD:DP:MM:MQS:NQSBQ:NQSMM:REnd:RStart:SC
	if(/^#CHROM/)
	{
		@sampleName=@F[-2,-1];
	}
	if(/^#/){ print "$_\n"; next; }
	if(!/SOMATIC/) { next; }
	my ($REnd,$RStart,$MM1,$MM2,$MQS1,$MQS2,$NQSBQ1,$NQSBQ2,$NQSMM1,$NQSMM2,$TAD1,$TAD2,$NAD1,$NAD2,$TSC,$NSC);
	my @K=split /:/,$F[-3];
	if($auto==0)
	{
		for(my $i=0;$i<@K;$i++)
		{
			if($K[$i] eq "GT")
			{
				my $_gt1=(split /:/,$F[-2])[$i];
				my $_gt2=(split /:/,$F[-1])[$i];
				if($_gt1 eq "0/0" && ($_gt2 eq "0/1" || $_gt2 eq "1/1"))
				{
					$opt_NIndex=-2;
					$opt_TIndex=-1;
					$auto=1;
					printf STDERR "$sampleName[-2],$sampleName[-1]\n";
					printf STDERR "NIndex:$opt_NIndex\tTIndex:$opt_TIndex\n";
				}elsif($_gt2 eq "0/0" && ($_gt1 eq "0/1" || $_gt1 eq "1/1"))
				{
					$opt_NIndex=-1;
					$opt_TIndex=-2;
					$auto=1;
					printf STDERR "$sampleName[-2],$sampleName[-1]\n";
					printf STDERR "NIndex:$opt_NIndex\tTIndex:$opt_TIndex\n";
				}else
				{
					printf STDERR "Can't distinguish normal/tumor\t$line\n";
					exit;
				}
			}
		}
	}
	my @VNormal=split /:/,$F[$opt_NIndex];
	my @VTumor=split /:/,$F[$opt_TIndex];
	for(my $i=0;$i<@K;$i++)
	{
		if($K[$i] eq "REnd") { $REnd=$VTumor[$i]; }
		if($K[$i] eq "RStart") { $RStart=$VTumor[$i]; }
		if($K[$i] eq "MM") { ($MM1,$MM2)=split ",",$VTumor[$i]; }
		if($K[$i] eq "MQS") { ($MQS1,$MQS2)=split ",",$VTumor[$i]; }
		if($K[$i] eq "NQSBQ") { ($NQSBQ1,$NQSBQ2)=split ",",$VTumor[$i]; }
		if($K[$i] eq "NQSMM") { ($NQSMM1,$NQSMM2)=split ",",$VTumor[$i]; }
		if($K[$i] eq "AD") { ($TAD1,$TAD2)=split ",",$VTumor[$i]; ($NAD1,$NAD2)=split ",",$VNormal[$i]; }
		if($K[$i] eq "SC") { $TSC=$VTumor[$i]; $NSC=$VNormal[$i]; }
	}
	my ($median_REnd,$mad_REnd)=split /,/,$REnd;
	my ($median_RStart,$mad_RStart)=split /,/,$RStart;
	
	if($median_RStart<$RStartTh || $median_REnd<$REndTh) { printf STDERR "read position\tobs:$median_RStart\tth:$RStartTh\t$_\n"; next; }
	if($MM1>$MMTh || $MM2>$MMTh) { printf STDERR "MM\tobs:$MM1,$MM2\tth:$MMTh\t$_\n"; next; }
	if($MQS2<$MQSTh) { printf STDERR "MQS2\tobs:$MQS2\tth:$MQSTh\t$_\n"; next; }
	if($NQSBQ2<$NQSBQTh) { printf STDERR "NQSBQ2\tobs:$NQSBQ2\tth:$NQSBQTh\t$_\n"; next; } 
	if($NQSMM2>$NQSMMTh) { printf STDERR "NQSMM2\tobs:$NQSMM2\tth:$NQSMMTh\t$_\n"; next; }
	### depth
	if($TAD1=~/\d+/ && $TAD2=~/\d+/ && $NAD1=~/\d+/ && $NAD2=~/\d+/)
	{
		if(($TAD1+$TAD2)<$TDepthTh) { printf STDERR "TDepth\tobs:$TAD1,$TAD2\tth:$TDepthTh\t$_\n"; next; }
		if(($NAD1+$NAD2)<$NDepthTh) { printf STDERR "NDepth\tobs:$NAD1,$NAD2\tth:$NDepthTh\t$_\n"; next; }
		if($NAD2/($NAD1+$NAD2)>$NFreqTh) { printf STDERR "NFreq\tobs:%4.2f\tth:$NFreqTh\t$_\n",$NAD2/($NAD1+$NAD2); next; }
		if($TAD2/($TAD1+$TAD2)<$TFreqTh) { printf STDERR "TFreq\tobs:%4.2f\tth:$TFreqTh\t$_\n",$TAD2/($TAD1+$TAD2); next; }
	}elsif($TSC=~/\d/ && $NSC=~/\d/)
	{
		my @_t=split /,/,$TSC;
		my $_d=0;
		foreach (@_t) { $_d+=$_; }
		if($_d < $TDepthTh) { printf STDERR "TDepth\tobs:$_d\tth:$TDepthTh\t$_\n"; next; }
		if(($_t[0]+$_t[1])/$_d<$TFreqTh) { printf STDERR "TFREQ\t%4.2f\tth:$TFreqTh\t$_\n",($_t[0]+$_t[1])/$_d; next; }

		@_t=split /,/,$NSC;
		$_d=0;
		foreach (@_t) { $_d+=$_; }
		if($_d < $NDepthTh) { printf STDERR "NDepth\tobs:$_d\tth:$NDepthTh\t$_\n"; next; }
		if(($_t[0]+$_t[1])/$_d>$NFreqTh) { printf STDERR "NFREQ\t%4.2f\tth:$NFreqTh\t$_\n",($_t[0]+$_t[1])/$_d; next; }

	}
	print "$_\n";
}



############################################################################
sub usage
{
	die `pod2text $0`;
}
