#!/usr/bin/perl
#============================================================================
# Name        		: varScanSomatic2annoVar.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Wed Apr 20 01:16:21 2011
# Last Modified By	: 
# Last Modified On	: Wed Apr 20 01:16:21 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl varScanSomatic2annoVar.pl [option] [infile]

	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }

while(<>)
{
	chomp;
	if(/^\s*$/) { next; }
	s/%//g;
	my @field=split /\t/;
	my ($chr,$pos,$ref,$alt,$gt_n,$gt_t)=@field[0,1,2,3,7,11];
	$field[6]=$field[6]/100;
	$field[10]=$field[10]/100;
	my $begin=$pos;
	my $end=$pos;
	if($alt=~/^\+(.+)/)
	{
		#ins
		$ref="-";
		$alt=$1;
	}elsif($alt=~/^\-(.+)/)
	{	#del
		$begin=$pos+1;
		$end=$pos+length($1);
		$ref=$1;
		$alt="-";
	}
	my $info=sprintf "TumorGT=%s;NormalReads=%s;TumorReads=%s;VariantP=%s;SomaticP=%s;TumorDP4=%s;NormalDP4=%s",$gt_t,join(",",@field[4,5,6]),join(",",@field[8,9,10]),$field[13],$field[14],join(",",@field[15..18]),join(",",@field[19..22]);
	#my @_info=();
	#for (my $i=0;$i< @info; $i++)
	#{
	#	push @_info,sprintf "%s=%s;",$info[$i],$field[$i+4];
	#}
	#if($gt_t=~/[ATCG]/) { $gt_t="hom"; }
	#else{ $gt_t="het"; }
	printf "$chr\t$begin\t$end\t$ref\t$alt\t$info\n";
	#chrom   position        ref     var     normal_reads1   normal_reads2   normal_var_freq normal_gt       tumor_reads1    tumor_reads2    tumor_var_freq  tumor_gt        somatic_status  variant_p_value somatic_p_value tumor_reads1_plus       tumor_reads1_minus      tumor_reads2_plus       tumor_reads2_minus
	#chr10   5993470 A       G       19      0       0%      A       11      3       21.43%  R       Somatic 1.0     0.0667155425219932      5       6       1       2
	#chr10   272610  G       -T      14      0       0%      G       7       2       22.22%  */-T    Somatic 1.0     0.14229249011857667     2       5       1       1
	#chr10   7302795 T       +TTC    10      0       0%      T       12      3       20%     */+TTC  Somatic 1.0     0.19782608695652215     6       6       1       2
}
############################################################################
sub usage
{
	die `pod2text $0`;
}
