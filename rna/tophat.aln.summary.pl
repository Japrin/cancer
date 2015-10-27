#!/usr/bin/perl
#============================================================================
# Name        		: tophat.aln.summary.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Wed Oct 15 11:24:48 2014
# Last Modified By	: 
# Last Modified On	: Wed Oct 15 11:24:48 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl tophat.aln.summary.pl [option] <infile>

    -s  sample [default: SAMPLE]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_s);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_s)) { $opt_s="SAMPLE"; }

my $totalRead=0;
my $mappedRead=0;
my $mappedPair=0;
my $mappingRate=0;
my $mappingConcordantRate=0;
my $concordantPair=0;
my $disconcordantPair=0;
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    
    if(/Input\s+:\s+(\d+)/){ $totalRead+=$1; };
    if(/Mapped\s+:\s+(\d+)/){ $mappedRead+=$1; };

    if(/Aligned pairs:\s+(\d+)/) { $mappedPair=$1; }
    if(/^(.+?) overall read mapping rate/) { $mappingRate=$1; }
    if(/^(.+?) concordant pair alignment rate/) { $mappingConcordantRate=$1; }
    if(/(\d+) (.+) are discordant alignments/) { $disconcordantPair=$1; }
}
$concordantPair=$mappedPair-$disconcordantPair;
print "#Sample\tTotalRead\tMappedRead\tMappedPair\tOverallMappingRate\tConcordantPair\tConcordantRate\n";
print "$opt_s\t$totalRead\t$mappedRead\t$mappedPair\t$mappingRate\t$concordantPair\t$mappingConcordantRate\n";
###################################################################




sub usage
{
    die `pod2text $0`;
}
