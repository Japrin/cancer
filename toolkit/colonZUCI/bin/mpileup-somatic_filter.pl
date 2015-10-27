#!/usr/bin/perl
#============================================================================
# Name        		: mpileup-somatic_filter.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Thu Sep 22 14:18:37 2011
# Last Modified By	: 
# Last Modified On	: Thu Sep 22 14:18:37 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl mpileup-somatic_filter.pl [option] <infile>

	--snv		filter those in dbSNP or in 1KG
	--CLR=<int>	CLR must >=<int>. default 15
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_snv,$opt_CLR);
GetOptions("h"	=>\$opt_h, "snv"=>\$opt_snv,"CLR"=>\$opt_CLR);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;
if(!defined($opt_CLR)) { $opt_CLR=15; }

if($infile =~ /\.gz$/) { open $in,"bgzip -cd $infile | " or die "Cann't open file $infile ($!)\n"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
while(<$in>)
{
	chomp;
	if(/^\s*$/ || /^#/) { printf "$_\n";next; }
	my @field=split /\t/;
	my ($chr,$pos,$rsID,$ref,$alt,$qual,$filter)=@field[0,1,2,3,4,5,6];
	if($opt_snv && $rsID ne ".") { next; }
	if($opt_snv && /1000G.*?=(.+?);/) { printf STDERR "in 1KG\tfreq:$1\t$_\n";next; }

	my ($Normal_ML_ALT_Dep,$Tumor_ML_ALT_Dep,$NormalMQ0,$TumorMQ0,$CLR,$normalDepth,$tumorDepth);
	($Normal_ML_ALT_Dep)= /Normal_ML_ALT_Dep=(.+?);/;
	($Tumor_ML_ALT_Dep)= /Tumor_ML_ALT_Dep=(.+?);/;
	($NormalMQ0)= /NormalMQ0=(.+?);/;
	($TumorMQ0)= /TumorMQ0=(.+?);/;
	($CLR)= /CLR=(.+?);/;
	my ($strandednessP)=/TumorStrandBias=(.+?);/;
	$normalDepth=(split /:/,$field[-2])[2];
	$tumorDepth=(split /:/,$field[-1])[2];
	if(!defined($CLR)) 
	{
		printf STDERR "NO CLR (probably no coverage in normal or tumor\t$_\n";
		next;
	}

	#chrX    128977251       .       G       T       11.3    .       Func=exonic;Gene=BCORL1;ExonicFunc=nonsynonymous SNV;AAChange=NM_021946:c.G2822T:p.G941V;SIFT=0;PolyPhen2=0.042;LJB_PhyloP=0.987207;LJB_MutationTaster=2.46E-4;LJB_LRT=0.956899;tfbs=Score=996(V$MZF1_01);ensGene=ENSG00000085185:exonic:nonsynonymous SNV:ENST00000218147:c.G2822T:p.G941V;DP=45;AF1=0.25;CLR=44;AC1=1;DP4=15,12,0,5;MQ=46;FQ=12.5;PV4=0.046,0.0029,0.068,0.04;ML_ALT=T;TumorRefDepth=7;Tumor_ML_ALT_Dep=5;NormalRefDepth=20;Normal_ML_ALT_Dep=0;NormalRefBQMean=24.00;NormalAltBQMean=0.00;TumorRefBQMean=24.29;TumorAltBQMean=19.40;NormalRefMQMean=52.09;NormalAltMQMean=0.00;TumorRefMQMean=37.00;TumorAltMQMean=37.00;NormalMQ0=0;TumorMQ0=0;NormalIndel=0;TumorIndel=0;DIST_MEAN=44.20;NormalStrandBias=1.0000;TumorStrandBias=0.4697;NormalBQBias=-1.0000;TumorBQBias=0.1579;NormalMQBias=-1.0000;TumorMQBias=0.8154;frequencyP=0.003933;       GT:PL:DP:SP:GQ  0/0:0,60,219:20:0:62    0/1:44,0,99:12:16:42
	
#	Normal_ML_ALT_Dep <= 1
#	|   Tumor_ML_ALT_Dep <= 3: No (6)
#	|   Tumor_ML_ALT_Dep > 3
#	|   |   NormalMQ0+TumorMQ0>=2: No (4/1)
#	|   |   else
#	|   |    |    CLR<18 : No (1)
#	|   |    |    CLR>=18: Yes (58/2)
#	Normal_ML_ALT_Dep > 1: No (7)
	#if(/4660017/) { print STDERR "4660017\tnormalDepth:$normalDepth\ttumorDepth:$tumorDepth\tstrandbias:$strandednessP\n"; }
	if($normalDepth<6) { next; }
	if($tumorDepth<6) { next; }
	if($strandednessP<0.01) { next; }

	if($Normal_ML_ALT_Dep>1) { next; }
	if($Tumor_ML_ALT_Dep<=3) { next; }
	if($NormalMQ0+$TumorMQ0>=2) { next; }
	if($CLR<$opt_CLR) { next; }
	printf "%s\n",$_;
}

############################################################################
sub usage
{
	die `pod2text $0`;
}
