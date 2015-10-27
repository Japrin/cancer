#!/usr/bin/env perl
#============================================================================
# Name        		: vcf.SampleReorder.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Tue Jan 13 08:43:10 2015
# Last Modified By	: 
# Last Modified On	: Tue Jan 13 08:43:10 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl vcf.SampleReorder.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;


if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
my @H=();
my $flag_h=0;
my $flag_f=0;
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^##/) { print "$line\n"; next; }
    my @F=split /\t/;
    if(/^#CHROM/) { @H=@F; next; }
    my @fmt=split /:/,$F[-3];
    my @sample1=split /:/,$F[-2];
    my @sample2=split /:/,$F[-1];
    my $iGT="";
    for(my $i=0;$i<@fmt;$i++)
    {
    	if($fmt[$i] eq "GT")
	{
		$iGT=$i;
		last;
	}
    }
    my $gt_sample1=$sample1[$iGT];
    my $gt_sample2=$sample2[$iGT];
    if(!$flag_f && ($gt_sample2 eq "0" || $gt_sample2 eq "0/0"))
    {
    	$flag_f=1;
    }
    if($flag_f==1)
    {
    	my $tmp=$F[-2];
	$F[-2]=$F[-1];
	$F[-1]=$tmp;
    }
    if(!$flag_h)
    {
    	if($flag_f)
	{
	    my $tmp=$H[-2];
	    $H[-2]=$H[-1];
	    $H[-1]=$tmp;
	}
	print join("\t",@H)."\n";
    	$flag_h=1;
    }
    print join("\t",@F)."\n";
}
###################################################################



sub usage
{
    die `pod2text $0`;
}
