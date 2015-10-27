#!/usr/bin/perl
#============================================================================
# Name        		: MentoCaloSpectrumDist.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Fri Sep  9 15:02:00 2011
# Last Modified By	: 
# Last Modified On	: Fri Sep  9 15:02:00 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl MentoCaloSpectrumDist.pl [option] <SiteDir> <spectrumFile>

	-i	simulation times [default 1000]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Dumper;
use Cwd qw(abs_path);
use Math::Random;

my ($in,$out);
my ($opt_h,$opt_i);
GetOptions("h"	=>\$opt_h,"i=i"=>\$opt_i);
if(@ARGV<2 || $opt_h) { usage(); }
my $siteDir=shift @ARGV;
my $spectrumFile=shift @ARGV;
$siteDir=abs_path $siteDir;
if(!defined($opt_i)) { $opt_i=1000; }

#A->T    0
#T->A    0
#A->C    6
#T->G    3
#A->G    6
#T->C    9
#C->A    11
#G->T    4
#C->TinCpG       20
#G->AinCpG       24
#C->TnotCpG      5
#G->AnotCpG      7
#C->G    0
#G->C    0

my %spectrum=();
readSpectrum(\%spectrum,$spectrumFile);

#A.bed.gz	8501657
#CinCpG.bed.gz	1119299
#CnotCpG.bed.gz	8129989
#GinCpG.bed.gz	1119275
#GnotCpG.bed.gz	8112472
#T.bed.gz	8538345
my %fileSize=();
readFileSize(\%fileSize,"$siteDir/data.summary.txt");
my %result=();
for(my $j=0;$j<$opt_i;$j++)
{
	$result{$j}={};
	while(randomOnce($result{$j},\%spectrum)!=1){ printf STDERR "retry in $j th simulation\n"; $result{$j}={}; };
}
foreach (sort {$a<=>$b} keys %result)
{
	my $p=$result{$_};
	foreach my $_f (keys %$p)
	{
		my $pp=$p->{$_f};
		foreach my $_l (keys %$pp)
		{
			printf "$_f\t$_l\t%s\t%s\t$_\n",$pp->{$_l}->[0],$pp->{$_l}->[1];
		}
	}
}

############################################################################
sub randomOnce
{
	my ($pOut,$pSpe)=@_;
	foreach (keys %$pSpe)
	{
		my ($ref,$alt)=(split //)[0,3];
		for(my $i=0;$i<$pSpe->{$_};$i++)
		{
			my $_file="";
			if($_ eq "C->TinCpG")
			{
				$_file="CinCpG.bed.gz";
			}elsif($_ eq "G->AinCpG")
			{
				$_file="GinCpG.bed.gz";
			}elsif($_ eq "C->TnotCpG")
			{
				$_file="CnotCpG.bed.gz";
			}elsif($_ eq "G->AnotCpG")
			{
				$_file="GnotCpG.bed.gz";
			}else
			{
				if($ref eq "C")
				{
					my $_c1=$fileSize{"CinCpG.bed.gz"};
					my $_c2=$fileSize{"CnotCpG.bed.gz"};
					my $_c=$_c1+$_c2;
					my $_d=int(rand($_c))+1;
					if($_d<=$_c1) { $_file="CinCpG.bed.gz"; }
					else { $_file="CnotCpG.bed.gz"; }
				}elsif($ref eq "G")
				{
					my $_c1=$fileSize{"GinCpG.bed.gz"};
					my $_c2=$fileSize{"GnotCpG.bed.gz"};
					my $_c=$_c1+$_c2;
					my $_d=int(rand($_c))+1;
					if($_d<=$_c1) { $_file="GinCpG.bed.gz"; }
					else { $_file="GnotCpG.bed.gz"; }
				}else
				{
					$_file="$ref.bed.gz";
				}
			}
			my $_totalSite=$fileSize{$_file};
			#my $_l=random_uniform_integer(1,1,$_totalSite);
			### test
			#my @_seed=random_get_seed();
			#printf STDERR "seed\t%s\t$_file\t$_l\n",join("\t",@_seed);
			### test
			my $_l=int(rand($_totalSite))+1;
			#printf STDERR "$_l\t$_totalSite\t$_file\n";
			if(!exists($pOut->{$_file})) { $pOut->{$_file}={}; }
			if(exists($pOut->{$_file}->{$_l})) { return 0; }
			else { $pOut->{$_file}->{$_l}=[$ref,$alt]; }
		}
	}
	return 1;
}
sub readFileSize
{
	my $in;
	my ($list,$infile)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/) { next; }
		my @field=split /\t/;
		my ($key,$c)=@field[0,1];
		$list->{$key}=$c;
	}
}
sub readSpectrum
{
	my $in;
	my ($list,$infile)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/) { next; }
		my @field=split /\t/;
		my ($key,$c)=@field[0,1];
		$list->{$key}=$c;
	}
}
sub usage
{
	die `pod2text $0`;
}
