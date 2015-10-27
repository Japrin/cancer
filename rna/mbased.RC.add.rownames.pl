#!/usr/bin/perl
#============================================================================
# Name        		: mbased.RC.add.rownames.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Dec 22 16:47:01 2014
# Last Modified By	: 
# Last Modified On	: Mon Dec 22 16:47:01 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl mbased.RC.add.rownames.pl [option] <infile>

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
#my $infile=shift @ARGV;

my $preGene="";
my $preI=0;
while(<>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    #chr1	2103506	C	G	PRKCZ	52	39	49	0
    #chr1	3703710	G	A	LRRC47	4	12	3	13
    #chr1	3713027	G	T	LRRC47	10	8	0	23
    my ($chr,$pos,$gene)=@F[0,1,4];
    if($gene ne $preGene)
    {
    	$preGene=$gene;
	$preI=0;
    }
    $preI++;
    print "$line\t${gene}_SNV$preI\n";
}
###################################################################




sub usage
{
    die `pod2text $0`;
}
