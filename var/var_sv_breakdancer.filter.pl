#!/usr/bin/env perl
#============================================================================
# Name        		: var_sv_breakdancer.filter.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Wed Dec  4 10:18:35 2013
# Last Modified By	: 
# Last Modified On	: Wed Dec  4 10:18:35 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl var_sv_breakdancer.filter.pl [option] <infile>

	-m		number of samples with sv [default: infer from header ]
	-n		number of supporting reads per sample [default 2]
	-g		gender [default "F", female]
	-a		only chromosome 1..22,X,Y
	-b		normalBam file [default OFF; use it when select somatic sv in normal/tumor pair]
	-minsize	minsize of deletion, insertion [default 50]
	-maxsize	maxsize of deletion, inversion [default 1e6]
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_n,$opt_m,$opt_g,$opt_a,$opt_minsize,$opt_maxsize,$opt_b);
GetOptions("h"	=>\$opt_h,"m=i"=>\$opt_m,"n=i"=>\$opt_n,"g=s"=>\$opt_g,"a"=>\$opt_a,"minsize=i"=>\$opt_minsize,"maxsize=i"=>\$opt_maxsize,"b=s"=>\$opt_b);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
#if(!defined($opt_m)) { $opt_m=1; }
if(!defined($opt_n)) { $opt_n=2; }
if(!defined($opt_g)) { $opt_g="F"; }
if(!defined($opt_minsize)) { $opt_minsize=50; }
if(!defined($opt_maxsize)) { $opt_maxsize=1e6; }

my %CHROMOSOME=('1'=>1,'2'=>1,'3'=>1,'4'=>1,'5'=>1,'6'=>1,'7'=>1,'8'=>1,'9'=>1,'10'=>1,'11'=>1,'12'=>1,'13'=>1,'14'=>1,'15'=>1,'16'=>1,'17'=>1,'18'=>1,'19'=>1,'20'=>1,'21'=>1,'22'=>1,'X'=>1,'Y'=>1);

my %para_m=();
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
	chomp;
	my $line=$_;

	#/PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH21.final.bam  mean:290.74     std:33.16       uppercutoff:405.92      lowercutoff
	if(/^#(.+?)\tmean:.+?\tstd:.+?\tuppercutoff/)
	{
		$para_m{$1}++;
	}
	if(/^#Chr1/)
	{
		if(!defined($opt_m))
		{
			$opt_m=keys %para_m;
			print STDERR "opt_m\t$opt_m\n";
		}
	}

	if(/^\s*$/ || /^#/) { print "$line\n";next; }
	my @F=split /\t/;
	#Chr1   Pos1    Orientation1    Chr2    Pos2    Orientation2    Type    Size    Score   num_Reads       num_Reads_lib   Allele_frequency        CH21.final.bam  CH22.final.bam  CH23.final.bam
	#1       10016   16+12-  1       10381   16+12-  ITX     -235    99      2       /PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH22.final.bam|2 -nan    NA      NA      NA
	#1       965794  11+14-  1       966028  1+3-    DEL     151     99      3       /PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH21.final.bam|1:/PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH22.final.bam|1:/PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH23.final.bam|1 0.89    0.25    0.19    0.25
	
	my ($chr1,$pos1,$chr2,$pos2)=@F[0,1,3,4];
	$chr1=~s/^chr//;
	$chr2=~s/^chr//;

	### a bug ? ###
	if($chr1 ne $chr2 && $F[6] ne "CTX")
	{
		print STDERR "WARNING: $chr1:$chr2:$F[6]\t$line\n";
		next;
	}
	
	### Filter by gender
	if($opt_g eq "F" || $opt_g=~/female/i)
	{
		if($chr1 eq "Y" || $chr2 eq "Y") { next; }
	}
	### Filter by chromosome
	if($opt_a)
	{
		if(!exists($CHROMOSOME{$chr1}) || !exists($CHROMOSOME{$chr2})) { next; }
	}
	
	my @num_reads_lib=split /:/,$F[10];
	
	my $num_sv_sample=scalar @num_reads_lib;
	
	my $isFiltered=0;
	if($opt_b)
	{
		## somatic mode: no sv supporting reads in normal, and supporting reads in tumor must >= opt_n
		foreach my $_l (@num_reads_lib)
		{
			my ($_s,$_n)=split /\|/,$_l;
			if($_s eq $opt_b) { $isFiltered=1; last; }
			elsif($_n<$opt_n) { $isFiltered=1; last; }
		}
		if($isFiltered) { next; }
	}else
	{
		## "shared" mode: number of samples with sv must >= opt_m && number of supporting reads per sample must >= opt_n
		if($num_sv_sample<$opt_m) { next; }
		foreach my $_l (@num_reads_lib)
		{
			my ($_s,$_n)=split /\|/,$_l;
			if($_n<$opt_n) { $isFiltered=1; last; }
		}
		if($isFiltered) { next; }
	}

	### Filter by size
	#-minsize	minsize of deletion, insertion [default 50]
	#-maxsize	maxsize of deletion, inversion [default 1e6]
	if($F[6] ne "ITX" && $F[6] ne "CTX")
	{
		my $_size= abs($pos1-$pos2);
		$_size= ($_size>abs($F[7])?$_size:abs($F[7]));
		if($_size<$opt_minsize || $_size>$opt_maxsize) { next; }
	}
	#if(($F[6] eq "DEL" || $F[6] eq "INS") && abs($_size)<$opt_minsize)
	#{
	#	next;
	#}
	#if(($F[6] eq "DEL" || $F[6] eq "INV") && abs($_size)>$opt_maxsize)
	#{
	#	next;
	#}

	### Pass all
	print "$line\n";

}
###################################################################




sub usage
{
	die `pod2text $0`;
}
