#!/usr/bin/perl
#============================================================================
# Name        		: somatic_cnv_varscan.gene.info.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Sat Dec 14 16:19:54 2013
# Last Modified By	: 
# Last Modified On	: Sat Dec 14 16:19:54 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl somatic_cnv_varscan.gene.info.pl [option] <infile>

	-s	sample
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_s);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s);
if(@ARGV<1 || $opt_h) { usage(); }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;

open $in,$infile or die "Cann't open file $infile ($!) \n";
print "Sample\tGene\tType\tseg_chr\tseg_begin\tseg_end\tseg_mean\tnum_segments\tnum_markers\tevent_size\tsize_class\tchrom_arm\tarm_fraction\tchrom_fraction\n";
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#2       232663578       233351301       ALPI,ALPP,ALPPL2,COPS7B,DIS3L2,ECEL1,NPPC;ECEL1P2,MIR1471       exonic;ncRNA_exonic     "2631748:L1HS(LINE)"    Score=0.942334;Name=chr2:233296814,chr2:233185679       deletion        .       seg_mean=-0.81325;num_segments=2;num_markers=198;event_size=687724;size_class=focal;chrom_arm=2q;arm_fraction=0.46%;chrom_fraction=0.28%
	if(/^#/ || /intergenic/){next;}
	my %g=();
	my @g=split /[,;]/,$F[3]; 
	my @info=split /[;=]/,$F[-1];
	my $info="";
	for(my $i=1;$i<=$#info;$i=$i+2)
	{
		#print STDERR "$i\t$#info\n";
		$info.="\t$info[$i]";
	} 
	foreach (@g)
	{
		$g{$_}=1;
	}
	foreach (keys %g)
	{
		print "$opt_s\t$_\t$F[7]\t$F[0]\t$F[1]\t$F[2]$info\n";
	}
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
