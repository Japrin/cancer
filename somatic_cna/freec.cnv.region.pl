#!/usr/bin/env perl
#============================================================================
# Name        		: freec.loh.region.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Fri Jan 10 00:26:32 2014
# Last Modified By	: 
# Last Modified On	: Fri Jan 10 00:26:32 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl freec.cnv.region.pl [option] <infile>
	
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_f);
GetOptions("h"	=>\$opt_h,"f=f"=>\$opt_f);
if(@ARGV<0 || $opt_h) { usage(); }
#my $infile=shift @ARGV;

my $preChr="";
my $preBeg="";
my $preEnd="";
my $preInfo="";
my $preCN="";
#open $in,$infile or die "Cann't open file $infile ($!) \n";
#my $h=<>;
while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#1       86005656        86008824        12      gain    AAAAAAAABBBB    -1      somatic 0
	#1       86008824        93890016        4       gain    AAAA    0.0426311       somatic 0
	if(@F<9) { next; }	
	my ($chr,$beg,$end,$CN)=@F[0,1,2,3];
	my $info=join("\t",@F[3,4,5],"NA",$F[7],"NA");

	if($chr eq $preChr && $beg==$preEnd && $CN==$preCN )
	{
		$preEnd=$end;
	}else
	{
		if($preChr ne "")
		{
			##output
			print "$preChr\t$preBeg\t$preEnd\t$preInfo\n";
		}
		$preChr=$chr;
		$preBeg=$beg;
		$preEnd=$end;
		$preInfo=$info;
		$preCN=$CN;
	}
}
if($preChr ne "")
{
	##output
	print "$preChr\t$preBeg\t$preEnd\t$preInfo\n";
}
###################################################################

sub usage
{
	die `pod2text $0`;
}
