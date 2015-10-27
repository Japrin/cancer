#!/usr/bin/perl
#============================================================================
# Name        		: mbased.addSNVDist.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Dec  8 17:11:59 2014
# Last Modified By	: 
# Last Modified On	: Mon Dec  8 17:11:59 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl mbased.addSNVDist.pl [option] <infile>
    -f	filter by "FILTER" [default OFF]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_a,$opt_b,$opt_f);
GetOptions("h"	=>\$opt_h,"a=i"=>\$opt_a,"b=f"=>\$opt_b,"f"=>\$opt_f);
if(@ARGV<0 || $opt_h) { usage(); }

my $lastLine="";
my $lastDist=300e6;
my $lastChr="";
my $lastPos="";

#if($infile=~/\.gz$/){ open $in,"bgzip -cd $infile |" or die "$!"; }
#else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
#while(<$in>)
while(<>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  466.466N
    #1       877715  rs6605066       C       G       266.77  PASS    Func=intronic;Gene=SAMD11;cpgIslandExt=(Name=CpG:_246);cytoband=1p36.33;1000g2012apr_all=0.89;snp138=rs6605066;ABHom=1.00;AC=2;AF=1.00;AN=2;DP=13;Dels=0.00;FS=0.000;HaplotypeScore=0.5784;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=20.52;VariantType=SNP   GT:AD:DP:GQ:PL  1/1:0,13:13:30:295,30,0
    #1       877782  rs79037098      C       G       404.77  PASS    Func=intronic;Gene=SAMD11;cpgIslandExt=(Name=CpG:_246);cytoband=1p36.33;1000g2012apr_all=0.07;snp138=rs79037098;ABHet=0.333;AC=1;AF=0.500;AN=2;BaseQRankSum=-1.694;DP=30;Dels=0.00;FS=0.000;HaplotypeScore=0.9789;MLEAC=1;MLEAF=0.500;MQ=60.00;MQ0=0;MQRankSum=0.242;QD=13.49;ReadPosRankSum=0.022;VariantType=SNP      GT:AD:DP:GQ:PL  0/1:10,20:30:99:433,0,231
    #
     
    if($opt_f && ( $F[6] ne "PASS" && $F[6] ne ".")) { next; }
    my ($chr,$pos)=@F[0,1];
    if($lastChr eq "")
    {
    	#first line
	$lastLine=$line;
	$lastChr=$chr;
	$lastPos=$pos;
	$lastDist=300e6;
    }else
    {
	($lastLine,$lastChr,$lastPos,$lastDist)=processOne($lastLine,$lastChr,$lastPos,$lastDist,$line,$chr,$pos);
    }
}
my @_F=split /\t/,$lastLine;
$_F[7].=";SNVDist=".$lastDist;
print join("\t",@_F)."\n";

###################################################################

sub processOne
{
	my ($lastLine,$lastChr,$lastPos,$lastDist,$line,$chr,$pos)=@_;
    	if($chr eq $lastChr)
	{
		my $newDist=abs($pos-$lastPos);
		if($newDist<$lastDist) { $lastDist=$newDist; }
		my @_F=split /\t/,$lastLine;
		$_F[7].=";SNVDist=".$lastDist;
		print join("\t",@_F)."\n";
		$lastLine=$line;
		$lastChr=$chr;
		$lastPos=$pos;
		$lastDist=$newDist;
	}else
	{
		my @_F=split /\t/,$lastLine;
		$_F[7].=";SNVDist=".$lastDist;
		print join("\t",@_F)."\n";
		$lastLine=$line;
		$lastChr=$chr;
		$lastPos=$pos;
		$lastDist=300e6;
		
	}
	return ($lastLine,$lastChr,$lastPos,$lastDist);

}

sub usage
{
    die `pod2text $0`;
}
