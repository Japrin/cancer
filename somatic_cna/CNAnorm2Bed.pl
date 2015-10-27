#!/usr/bin/perl
#============================================================================
# Name        		: CNAnorm2Bed.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Sat May  4 18:43:45 2013
# Last Modified By	: 
# Last Modified On	: Sat May  4 18:43:45 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl CNAnorm2Bed.pl [option] [infile]

	-w	window size [required]
	-show	[default ploidy, others: ratio]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_w,$opt_show);
GetOptions("h"	=>\$opt_h,"w=i"=>\$opt_w,"show=s"=>\$opt_show);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_w)) { usage(); }
if(!defined($opt_show)) { $opt_show="ploidy"; }
if($opt_show ne "ploidy" && $opt_show ne "ratio") { usage(); }

my $header=<>;
my $segCN="";
my $segChr="";
my $segBeg="";
my $segEnd="";
my $segRatio="";
while(<>)
{
	chomp;
	s/"//g;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#"Chr"   "Pos"   "Ratio" "Ratio.n"       "Ratio.s.n"     "SegMean"       "SegMean.n"
	#"chr1"  1       2       4.44320010183847        4.33803085153195        1.82867665301389        3.92763060495256
	#"chr1"  26001   2.5     5.94786804685502        4.27949787257732        1.82867665301389        3.92763060495256
	#"chr1"  52001   3       7.45253599187158        4.22099491432524        1.82867665301389        3.92763060495256
	if($F[-1] eq $segCN && $F[0] eq $segChr)
	{
		$segEnd=$F[1]+$opt_w-1;
	}else
	{
		if($segChr ne "")
		{
			#chr     beg     end     length  gain/loss       ratio   T_reads N_reads other1 other2
			my $_flag="";
			if($opt_show eq "ploidy") { $_flag=$segCN eq "NA"?"NA":($segCN>2?"gain":"loss"); }
			else{ $_flag=$segCN eq "NA"?"NA":($segCN>1?"gain":"loss"); }
			printf "$segChr\t%d\t$segEnd\t%d\t%s\t%s\tNA\tNA\t$segCN\tNA\n",$segBeg-1,$segEnd-$segBeg+1,$_flag,$segRatio;
		}
		$segCN=$F[-1];
		$segRatio=$F[-2];
		$segChr=$F[0];
		$segBeg=$F[1];
		$segEnd=$segBeg+$opt_w-1;
	}
}
if($segChr ne "")
{
	my $_flag="";
	if($opt_show eq "ploidy") { $_flag=$segCN eq "NA"?"NA":($segCN>2?"gain":"loss"); }
	else{ $_flag=$segCN eq "NA"?"NA":($segCN>1?"gain":"loss"); }
	printf "$segChr\t%d\t$segEnd\t%d\t%s\t%s\tNA\tNA\t$segCN\tNA\n",$segBeg-1,$segEnd-$segBeg+1,$_flag,$segRatio;
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
