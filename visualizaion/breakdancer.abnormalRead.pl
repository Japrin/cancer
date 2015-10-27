#!/usr/bin/perl
#============================================================================
# Name        		: breakdancer.abnormalRead.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Tue Sep 24 19:13:47 2013
# Last Modified By	: 
# Last Modified On	: Tue Sep 24 19:13:47 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl breakdancer.abnormalRead.pl [option] <configure> <inbam> 

	-r	region used for samtools view [required]
	-q	minimum mapping quality [default 30]
	-d	only Deletion support read pairs [default OFF]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_r,$opt_q,$opt_d);
GetOptions("h"	=>\$opt_h,"r=s"=>\$opt_r,"q=i"=>\$opt_q,"d"=>\$opt_d);
if(@ARGV<2 || $opt_h) { usage(); }
if(!defined($opt_r)) { usage(); }
if(!defined($opt_q)) { $opt_q=30; }
my $confFile=shift @ARGV;
my $infile=shift @ARGV;

my %conf=();
readList(\%conf,$confFile);
my %lib=();
readBamHeader(\%lib,$infile);

my %reads=();
$opt_r=~s/^chr//;
open $in,"samtools view -X -q $opt_q $infile $opt_r |" or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#HWI-ST1276:121:C28LNACXX:5:1108:19112:33413     163     1       9996    23      74M26S  =       10022   112     TTCCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCTTAACCCTAACC    >?>@A@A;@A>=A@=><<A@>?=AA@@@@AA>?=ABB@>@A@@@@7:.;,>==BB?@@>B?;><@C<:<@BBBC##########################    X0:i:1  X1:i:1  XA:Z:4,+10114,100M,3;   XC:i:74 MD:Z:1C1G70     PG:Z:MarkDuplicates     RG:Z:CS-3       XG:i:0  AM:i:0  NM:i:2  SM:i:23 XM:i:2  XN:i:5  XO:i:0  MQ:i:7  OQ:Z:B@@FFDD;FFAAFGDEE>GCFGEGGGIIGGCBD@FGIG@FH?BFB28(9)888=C@H@9D@7;6AA77;;?>AC##########################       XT:A:U
	my $_id;
	for(my $i=11;$i<@F;$i++)
	{
		if($F[$i]=~/^RG:Z:(.+)$/) { $_id=$1; last; }
	}
	my $_lib=$lib{$_id};
	#print "$_id\t$_lib\n";
	my ($_upper,$_lower)=@{$conf{$_lib}};
	my $isize=abs($F[8]);
	my $rid=$F[0];
	my $isAB=0;
	if($opt_d && $_upper<$isize && $F[2] eq $F[7])
	{
		$isAB=1;
	}elsif(!$opt_d && $isize<$_lower || $_upper<$isize)
	{
		$isAB=1;
	}
	if($isAB)
	{
		printf STDERR "$line\n";
		if(exists($reads{$rid}))
		{
			#
		}else
		{
			printf "$rid\t$F[2]\t$F[3]\t%s\t$F[7]\t$isize\n",$F[6] eq "="?$F[2]:$F[6];
			$reads{$rid}=$isize;
		}
	}
}
printf STDERR "Total read pairs:\t%d\n",scalar (keys %reads);
###################################################################

sub readBamHeader
{
	my $in;
	my ($pList,$infile)=@_;
	open $in,"samtools view -H $infile |" or die "Cann't open file $infile ($!) \n";
	my $h=<$in>;
	while(<$in>)
	{
		chomp;
		my $line=$_;
		my @F=split /\t/;
		if(/^\@RG/)
		{
			#@RG     ID:CS-3.3       PL:illumina     PU:NHD0677_L2   LB:NHD0677      SM:CS-3 CN:novogene
			my ($id,$library);
			for(my $i=1;$i<@F;$i++)
			{
				if($F[$i]=~/^ID:(.+)$/) { $id=$1; }
				if($F[$i]=~/^LB:(.+)$/) { $library=$1; }
			}
			$pList->{$id}=$library;
		}
	}
}

sub readList
{
	my $in;
	my ($pList,$infile)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	my $h=<$in>;
	while(<$in>)
	{
		chomp;
		my $line=$_;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
		#library uppercutoff     lowercutoff
		#NHD0677 358.21  175.01
		#NHD0695 596.24  0
		$pList->{$F[0]}=[$F[1],$F[2]];
	}
}

sub usage
{
	die `pod2text $0`;
}
