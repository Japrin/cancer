#!/usr/bin/perl
#============================================================================
# Name        		: stat.varScan.cnv.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Sun Jun  9 15:55:17 2013
# Last Modified By	: 
# Last Modified On	: Sun Jun  9 15:55:17 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl stat.varScan.cnv.pl [option] [infile]

	-c	not sex chromosome
	-p	prefix [defualt ./]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Dumper;

my ($in,$out);
my ($opt_h,$opt_c,$opt_p);
GetOptions("h"	=>\$opt_h,"c"=>\$opt_c,"p=s"=>\$opt_p);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_p)) { $opt_p="./"; }
#my $outfile=shift @ARGV;
#my $infile=shift @ARGV;

#open $in,$infile or die "Cann't open file $infile ($!) \n";
my $bin_len_all=0;
my $bin_len_gain=0;
my $bin_len_loss=0;
my $bin_count_all=0;
my $bin_count_gain=0;
my $bin_count_loss=0;

my $seg_len_all=0;
my $seg_len_gain=0;
my $seg_len_loss=0;
my %seg_count_all=();
my %seg_count_gain=();
my %seg_count_loss=();

my %gene_all=();
my %gene_gain=();
my %gene_loss=();

while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	if($opt_c && $F[0]=~/X|Y/) { next; }
	$bin_len_all+=$F[3];
	#$seg_len_all+=$F[12]-$F[11];
	$seg_count_all{$F[16]}=$F[12]-$F[11];
	$bin_count_all++;
	#PRDM16:NM_022114:Exon15;PRDM16:NM_199454:Exon15
	my @gg=();
	if($F[20] ne ".")
	{
		my @g=split /;/,$F[20];
		foreach (@g)
		{
			my $_g=(split /:/)[0];
			push @gg,$_g;
			$gene_all{$_g}=1;
		}
	}
	if($F[15] eq "gain")
	{
		$bin_len_gain+=$F[3];
		#$seg_len_gain+=$F[12]-$F[11];
		$seg_count_gain{$F[16]}=$F[12]-$F[11];
		$bin_count_gain++;
		foreach (@gg) { $gene_gain{$_}=1; }
	}elsif($F[15] eq "loss")
	{
		$bin_len_loss+=$F[3];
		#$seg_len_loss+=$F[12]-$F[11];
		$seg_count_loss{$F[16]}=$F[12]-$F[11];
		$bin_count_loss++;
		foreach (@gg) { $gene_loss{$_}=1; }
	}
}
foreach (keys %seg_count_all) { $seg_len_all+=$seg_count_all{$_}; }
foreach (keys %seg_count_gain) { $seg_len_gain+=$seg_count_gain{$_}; }
foreach (keys %seg_count_loss) { $seg_len_loss+=$seg_count_loss{$_}; }

open $out,">","$opt_p.summary.txt" or die "$!";
printf $out "bin count\t$bin_count_all\n";
printf $out "bin count (gain)\t$bin_count_gain\n";
printf $out "bin count (loss)\t$bin_count_loss\n";
printf $out "bin length\t$bin_len_all bp\n";
printf $out "bin length (gain)\t$bin_len_gain bp\n";
printf $out "bin length (loss)\t$bin_len_loss bp\n";
printf $out "segment count\t%d\n",scalar (keys %seg_count_all);
printf $out "segment count (loss)\t%d\n",scalar (keys %seg_count_loss);
printf $out "segment count (gain)\t%d\n",scalar (keys %seg_count_gain);
printf $out "segment length\t$seg_len_all bp\n";
printf $out "segment length (gain)\t$seg_len_gain bp\n";
printf $out "segment length (loss)\t$seg_len_loss bp\n";
printf $out "gene\t%d\n",scalar (keys %gene_all);
printf $out "gene (gain)\t%d\n",scalar (keys %gene_gain);
printf $out "gene (loss)\t%d\n",scalar (keys %gene_loss);

open $out,">","$opt_p.gain.gene.list" or die "$!";
foreach (keys %gene_gain) { print $out "$_\n"; }
open $out,">","$opt_p.loss.gene.list" or die "$!";
foreach (keys %gene_loss) { print $out "$_\n"; }

###################################################################




sub usage
{
	die `pod2text $0`;
}
