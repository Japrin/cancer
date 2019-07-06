#!/usr/bin/perl
#============================================================================
# Name        		: XXX2fa.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Oct 13 10:41:57 2014
# Last Modified By	: 
# Last Modified On	: Mon Oct 13 10:41:57 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl XXX2fa.pl [option] <infile>
    -o  output prefix [default: ./out ]
    -l  pep length [default: 9]
    -m  max output sequence length; NULL to disable [default 20000]
    -s  output sequence name's prefix  [default "sample"]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path qw(make_path remove_tree);

my ($in,$out,$out1,$out2,$out3,$out4);
my ($opt_h,$opt_e,$opt_o,$opt_l,$opt_m,$opt_s,$opt_i);
GetOptions("h"	=>\$opt_h,"e=s"=>\$opt_e,"o=s"=>\$opt_o,"l=i"=>\$opt_l,"m=i"=>\$opt_m,"s=s"=>\$opt_s,"i=s"=>\$opt_i);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_o)) { $opt_o="./out"; }
if(!defined($opt_l)) { $opt_l=9; }
if(!defined($opt_m)) { $opt_m=20000; }
if(!defined($opt_s)) { $opt_s="sample"; }

my ($outDir)=$opt_o=~/(.+)\//;
make_path($outDir);

open $in,$infile or die "Cann't open file $infile ($!) \n";
open $out1,">","$opt_o.fa" or die "$!";
open $out2,">","$opt_o.bed" or die "$!";

my $gID=1;
## gBeg, gEnd: bed format coordinate
my $gBeg=0;
my $gEnd=0;
my $gSeq1="";
my $gSeq2="";

while(<$in>)
{
    chomp;
    my $ID=$_;
	$ID=~s/^>//;
	my $seq=<$in>;
	chomp $seq;
	my $pipetideLen=length($seq);

	my $o_name=$ID;
	#my $o_beg=-1;
	#my $o_end=-1;

	if($opt_m ne "NULL" && $gEnd+$pipetideLen>$opt_m)
	{
		#print $out1 ">${opt_s}_${gID}\n$gSeq1\n";
		print $out1 ">${gID}\n$gSeq1\n";
		$gID++;
		$gSeq1=$seq;
		$gBeg=0;
	}else
	{
		$gBeg=$gEnd;
		$gSeq1.=$seq;
	}
	$gEnd=$gBeg+$pipetideLen;
	print $out2 "$gID\t".$gBeg."\t".$gEnd."\t$o_name\n";
}
if($gSeq1 ne "")
{
	print $out1 ">${gID}\n$gSeq1\n";
}
###################################################################

sub addPadding
{
	my ($seq,$mutPos,$peptideLen)=@_;
	###  peptideLen=9
	###  length(seq)=14
	###  mutPos=8
	###  -------X------
	### . n1=1
	###            x=6
	###                ..  n2=2  
	my $n1=$peptideLen-$mutPos;
	my $x=length($seq)-$mutPos;
	my $n2=$peptideLen-1-$x;
	return ("A"x$n1.$seq."A"x$n2,$n1,$n2);
}

sub usage
{
    die `pod2text $0`;
}
