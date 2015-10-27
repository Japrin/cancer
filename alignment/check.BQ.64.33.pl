#!/usr/bin/perl
#============================================================================
# Name        		: check.BQ.64.33.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Dec 29 18:28:58 2014
# Last Modified By	: 
# Last Modified On	: Mon Dec 29 18:28:58 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl check.BQ.64.33.pl [option] <infile>

    -o	verbose output prefix, such as "output"; if not specified, no output
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out,$out_1,$out_2);
my ($opt_h,$opt_o);
GetOptions("h"	=>\$opt_h,"o=s"=>\$opt_o);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
#if(!defined($opt_o)) { $opt_o="output"; }

my $is_64=1;
my $is_33=1;

my $i=0;
my $j=0;
if($infile=~/\.gz$/) { open $in,"gzip -cd $infile |" or die "$!"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }

if($opt_o)
{
    open $out_1,">","$opt_o.verbose.64.log" or die "$!";
    open $out_2,">","$opt_o.verbose.33.log" or die "$!";
}

while(<$in>)
{
    chomp;
    my $line=$_;
    #my @F=split /\t/;
    $i++;
    if($i%4==0)
    {
    	$j++;
    	my @a=split //,$_;
	my @b=map { ord($_)-64 } @a;
	my @c=map { ord($_)-33 } @a;
	foreach (@b)
	{
		if($_ > 42 || $_ <0) { $is_64=0; }
	
	}
	foreach (@c)
	{
		if($_ > 42 || $_ <0) { $is_33=0; }

	}
	
	if($opt_o)
	{
	    print $out_1 scalar @b."\t".join("\t",@b)."\n"; 
	    print $out_2 scalar @c."\t".join("\t",@c)."\n"; 
	}
	if($j>100) { last; }
    }
}

if($is_33 && !$is_64) { print "33\n"; }
elsif(!$is_33 && $is_64) { print "64\n"; }
else { print "NA\n"; }
###################################################################




sub usage
{
    die `pod2text $0`;
}
