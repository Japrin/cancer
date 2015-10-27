#!/usr/bin/env perl
#============================================================================
# Name        		: gistic.amp-del_genesToTable.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Sat Jan 17 21:10:38 2015
# Last Modified By	: 
# Last Modified On	: Sat Jan 17 21:10:38 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl gistic.amp-del_genesToTable.pl [option] <infile>

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
my @cytobandF=();
my @QValueF=();
my @resQF=();
my @widePickF=();

print "Gene Symbol\tcytoband\tq value\tresidual q value\twide peak boundaries\n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    s/\s\+$//;
    my @F=split /\t/;
    #cytoband        11q13.4 3q26.33 20q11.22
    #q value 2.2324e-05      0.063632        0.11117
    #residual q value        2.2324e-05      0.063632        0.11117
    #wide peak boundaries    chr11:68855231-73357576 chr3:158178735-183994731        chr20:33896835-35294588
    #genes in wide peak      hsa-mir-139     hsa-mir-1224    hsa-mir-1289-1
    #        hsa-mir-3165    hsa-mir-569     EPB41L1
    #        hsa-mir-548k    hsa-mir-551b    SPAG4
    if(/^cytoband/){ @cytobandF=@F; next; } 
    elsif(/^q value/){ @QValueF=@F; next; }
    elsif(/^residual q/) { @resQF=@F; next; }
    elsif(/^wide peak boundaries/) { @widePickF=@F; next; }
    for(my $i=1;$i<@F;$i++)
    {
    	if($F[$i] ne "")
	{
    		print "$F[$i]\t$cytobandF[$i]\t$QValueF[$i]\t$resQF[$i]\t$widePickF[$i]\n";
	}
    }
}
###################################################################




sub usage
{
    die `pod2text $0`;
}
