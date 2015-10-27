#!/usr/bin/perl
#============================================================================
# Name        		: var_sv_meerkat.format.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Sat Oct 11 15:12:52 2014
# Last Modified By	: 
# Last Modified On	: Sat Oct 11 15:12:52 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl var_sv_meerkat.format.pl [option] <infile>

    -p  id prefix [default "SVID"]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_p);
GetOptions("h"	=>\$opt_h,"p=s"=>\$opt_p);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_p)) { $opt_p="SVID"; }

my $id=0;
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    my $eventType=$F[0];
    $id++;
    if($eventType eq "del" || $eventType eq "tandem_dup" || $eventType=~/invers*/)
    {
        print "$F[5]\t$F[6]\t$F[7]\t0\t0\t$opt_p.$id\t$line\n";
    }elsif($eventType eq "del_invers" || $eventType =~/ins[so][ud]/)
    {
        print "$F[5]\t$F[6]\t$F[7]\t0\t0\t$opt_p.$id\t$line\n";
        print "$F[9]\t$F[10]\t$F[11]\t0\t0\t$opt_p.$id\t$line\n";
    }elsif($eventType =~/ins[so]$/)
    {
        print "$F[5]\t$F[6]\t$F[7]\t0\t0\t$opt_p.$id\t$line\n";
        print "$F[9]\t$F[10]\t$F[11]\t0\t0\t$opt_p.$id\t$line\n";
    }elsif($eventType eq "transl_inter")
    {
        print "$F[5]\t$F[6]\t$F[6]\t0\t0\t$opt_p.$id\t$line\n";
        print "$F[8]\t$F[9]\t$F[9]\t0\t0\t$opt_p.$id\t$line\n";
    }elsif($eventType eq "del_ins")
    {
        print "$F[5]\t$F[6]\t$F[7]\t0\t0\t$opt_p.$id\t$line\n";
    }
    else
    {
        print"ERROR\t$line\n";
    }
}
###################################################################




sub usage
{
    die `pod2text $0`;
}
