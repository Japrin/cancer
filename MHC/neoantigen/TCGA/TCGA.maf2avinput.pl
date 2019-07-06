#!/usr/bin/env perl
#============================================================================
# Name        		: TCGA.maf2avinput.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Tue Apr  7 21:28:00 2015
# Last Modified By	: 
# Last Modified On	: Tue Apr  7 21:28:00 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl TCGA.maf2avinput.pl [option] <infile>

    -i  columns for chr,pos,ref,alt, "," seperated and 0-based. [default 4,5,11,12 ]
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
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
}
###################################################################




sub usage
{
    die `pod2text $0`;
}
