#!/usr/bin/perl
#============================================================================
# Name        		: varScan-somatic_filter.pl
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

	perl varScan-somatic_filter.pl [option] <infile>

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
	if($opt_snv && /1000G.*?=(.+?);/) { print STDERR "in 1KG\tfreq:$1\t$_\n";next; }

#chr2    32494045        .       G       T       .       .       Func=exonic;Gene=BIRC6;ExonicFunc=stopgain SNV;AAChange=NM_016252:c.G2182T:p.E728X;Conserved=669(lod=678);SIFT=0;ensGene=ENSG00000115760:exonic:stopgain SNV:ENST00000261359:c.G2098T:p.E700X;normal_reads1=57;normal_reads2=0;normal_var_freq=0%;normal_gt=G;tumor_reads1=51;tumor_reads2=20;tumor_var_freq=28.17%;tumor_gt=K;somatic_status=Somatic;variant_p_value=1.0;somatic_p_value=1.8834613129936157E-6;tumor_reads1_plus=27;tumor_reads1_minus=24;tumor_reads2_plus=7;validation=Yes;
	
	#tumor_reads2 <= 3: No (12.0)
	#tumor_reads2 > 3
	#|   normal_var_freq <= 0.05
	#|   |   : No (/)
	#|   |   : Yes (/)
	#|   normal_var_freq > 0.05: No (3.0/1.0)
	my ($tumor_reads1,$tumor_reads2,$normal_reads1,$normal_reads2,$normal_var_freq,$tumor_var_freq);
	($tumor_reads1)=/tumor_reads1=(.+?);/;
	($tumor_reads2)=/tumor_reads2=(.+?);/;
	($normal_reads1)=/normal_reads1=(.+?);/;
	($normal_reads2)=/normal_reads2=(.+?);/;
	($tumor_var_freq)=/tumor_var_freq=(.+?);/;
	($normal_var_freq)=/normal_var_freq=(.+?);/;
	if($tumor_var_freq =~ /(.+)\%$/) { $tumor_var_freq=$1/100.0; }
	if($normal_var_freq =~ /(.+)\%$/) { $normal_var_freq=$1/100.0; }
	
	if($tumor_reads2<=3) { next; }
	#if($tumor_reads1+$tumor_reads2<6) { next; }
	#if($normal_reads1+$normal_reads2<6) { next; }
	if($normal_var_freq>0.05) { next; }
	#if($tumor_var_freq<0.2) { next; }

	printf "%s\n",$_;
}

############################################################################
sub usage
{
	die `pod2text $0`;
}
