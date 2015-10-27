#!/usr/bin/perl
#============================================================================
# Name        		: JointSNVMix-somatic_filter.pl
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

	perl JointSNVMix-somatic_filter.pl [option] <infile>

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
	#chr1    10308823        .       A       C       .       .       Func=exonic;Gene=KIF1B;ExonicFunc=nonsynonymous SNV;AAChange=NM_015074:c.A2605C:p.T869P;Conserved=672(lod=701);SIFT=0.27;PolyPhen2=0.949;LJB_PhyloP=0.998893;LJB_MutationTaster=0.943235;LJB_LRT=1.0;tfbs=Score=968(V$NFAT_Q6);ensGene=ENSG00000054523:exonic:nonsynonymous SNV:ENST00000263934:c.A2605C:p.T869P;somaticProbability=0.990366728835;NormalNonRef=C;TumorNonRef=C;NormalRefCounts=0;NormalNonRefCounts=0;TumorRefCounts=9;TumorNonRefCounts=4;p_AA_AA=0.00958679301502;p_AA_AB=0.990366728835;p_AA_BB=0.0;p_AB_AA=4.81280856702e-12;p_AB_AB=4.64781454951e-05;p_AB_BB=0.0;p_BB_AA=0.0;p_BB_AB=0.0;p_BB_BB=0.0;validation=NA;
	my ($somaticProbability,$NormalRefCounts,$NormalNonRefCounts,$TumorRefCounts,$TumorNonRefCounts,$NormalMQ0,$TumorMQ0);
	($somaticProbability)=/somaticProbability=(.+?);/;
	($NormalRefCounts)=/NormalRefCounts=(.+?);/;
	($NormalNonRefCounts)=/NormalNonRefCounts=(.+?);/;
	($TumorRefCounts)=/TumorRefCounts=(.+?);/;
	($TumorNonRefCounts)=/TumorNonRefCounts=(.+?);/;
	($NormalMQ0)=/NormalMQ0=(.+?);/;
	($TumorMQ0)=/TumorMQ0=(.+?);/;
	
	my ($strandednessP)=/TumorStrandBias=(.+?);/;
	if($strandednessP<0.05) { next; }

	if($NormalRefCounts+$NormalNonRefCounts<6) { next; }
	if($TumorRefCounts+$TumorNonRefCounts<6) { next; }
	if($somaticProbability<0.95) { next; }
	
	if($NormalNonRefCounts>1) { next; }
	if($TumorNonRefCounts<=3) { next; }
	if($NormalMQ0+$TumorMQ0>=2) { next; }
	printf "%s\n",$_;
}

############################################################################
sub usage
{
	die `pod2text $0`;
}
