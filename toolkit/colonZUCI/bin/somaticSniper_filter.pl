#!/usr/bin/perl
#============================================================================
# Name        		: somaticSniper_filter.pl
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

	perl somaticSniper_filter.pl [option] <infile>

	--snv	filter those in dbSNP or in 1KG
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_snv);
GetOptions("h"	=>\$opt_h, "snv"=>\$opt_snv);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;

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

	my ($Normal_ML_ALT_Dep,$Tumor_ML_ALT_Dep,$NormalMQ0,$TumorMQ0,$somaticScore,$DepthTumor,$DepthNormal);
	($Normal_ML_ALT_Dep)= /Normal_ML_ALT_Dep=(.+?);/;
	($Tumor_ML_ALT_Dep)= /Tumor_ML_ALT_Dep=(.+?);/;
	($NormalMQ0)= /NormalMQ0=(.+?);/;
	($TumorMQ0)= /TumorMQ0=(.+?);/;
	($somaticScore)= /somaticScore=(.+?);/;
	($DepthTumor)=/DepthTumor=(.+?);/;
	($DepthNormal)=/DepthNormal=(.+?);/;

#chr2    32494045        .       G       T       .       .       Func=exonic;Gene=BIRC6;ExonicFunc=stopgain SNV;AAChange=NM_016252:c.G2182T:p.E728X;Conserved=669(lod=678);SIFT=0;ensGene=ENSG00000115760:exonic:stopgain SNV:ENST00000261359:c.G2098T:p.E700X;somaticScore=109;tumorConsensusQuality=112;tumorSNVQuality=112;tumorRMSMappingQuality=37;DepthTumor=72;DepthNormal=63;validation=Yes;ML_ALT=T;TumorRefDepth=51;Tumor_ML_ALT_Dep=20;NormalRefDepth=59;Normal_ML_ALT_Dep=0;NormalRefBQMean=22.97;NormalAltBQMean=0.00;TumorRefBQMean=23.18;TumorAltBQMean=25.10;NormalRefMQMean=45.19;NormalAltMQMean=0.00;TumorRefMQMean=37.00;TumorAltMQMean=37.00;NormalMQ0=0;TumorMQ0=0;NormalIndel=0;TumorIndel=0;DIST_MEAN=45.00;NormalStrandBias=1.0000;TumorStrandBias=0.4308;NormalBQBias=-1.0000;TumorBQBias=0.0012;NormalMQBias=-1.0000;TumorMQBias=0.5765;frequencyP=0.000002;
#	Normal_ML_ALT_Dep <= 1
#	|   Tumor_ML_ALT_Dep <= 3: No (6)
#	|   Tumor_ML_ALT_Dep > 3
#	|   |   NormalMQ0+TumorMQ0>=2: No (4/1)
#	|   |   else
#	|   |    |    CLR<18 : No (1)
#	|   |    |    CLR>=18: Yes (58/2)
#	Normal_ML_ALT_Dep > 1: No (7)
#	my ($strTumorStrand)=/TumorStrand=(.+?);/;
#	my @aryTumorStrand=split /,/,$strTumorStrand;
#	if(($aryTumorStrand[0]>0 || $aryTumorStrand[1]>0) && ($aryTumorStrand[2]>0 || $aryTumorStrand[3]>0))
#	{
#		my $strandedness1=$aryTumorStrand[0]/($aryTumorStrand[0]+$aryTumorStrand[1]);
#		my $strandedness2=$aryTumorStrand[2]/($aryTumorStrand[2]+$aryTumorStrand[3]);
#		my $strandednessDiff=abs($strandedness1-$strandedness2);
#		if($strandednessDiff > 0.10 && ($strandedness2 < 0.10 || $strandedness2 > 0.90)) { next; }
#	}
	my ($strandednessP)=/TumorStrandBias=(.+?);/;
	if($strandednessP<0.05) { next; }

	if($Normal_ML_ALT_Dep>1) { next; }
	if($Tumor_ML_ALT_Dep<=3) { next; }
	if($NormalMQ0+$TumorMQ0>=2) { next; }
	if($DepthTumor<6) { next; }
	if($DepthNormal<6) { next; }
	#if($somaticScore<13) { next; }
	printf "%s\n",$_;
}

############################################################################
sub usage
{
	die `pod2text $0`;
}
