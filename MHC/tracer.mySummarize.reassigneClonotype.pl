#!/usr/bin/env perl
#============================================================================
# Name        		: tracer.mySummarize.reassigneClonotype.pl
# Author      		: zhenglt@gmail.com
# Version     		: v1.00
# Created On  		: Thu Dec 10 09:10:04 2015
# Last Modified By	: 
# Last Modified On	: Thu Dec 10 09:10:04 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl tracer.mySummarize.reassigneClonotype.pl [option] <infile> <outfile>

    -s  sampleID [default: SAMPLE]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_s);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s);
if(@ARGV<2 || $opt_h) { usage(); }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
my $outfile=shift @ARGV;

run($infile,$outfile);

###################################################################



sub run
{
    my ($in,$out);
    my ($infile,$outfile)=@_;
    my %clonotype=();
    my %clonotypeCount=();
    my @inputData=();
    my %cellClonotype=();
    open $in,$infile or die "Cann't open file $infile ($!) \n";
    open $out,">",$outfile or die "Cann't open file $outfile ($!) \n";
    my $header=<$in>;
    chomp $header;
    print $out "$header\tC_strict\tC_share\tC_strict_b\tC_share_b\n";
    my $cid=0;
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
        my $cellID=$F[0];
        my @idAlphaBeta4=@F[3,11,19,27];
        if($clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"}
            || $clonotype{"$idAlphaBeta4[1]:$idAlphaBeta4[0]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"}
            || $clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[3]:$idAlphaBeta4[2]"}
            || $clonotype{"$idAlphaBeta4[1]:$idAlphaBeta4[0]:$idAlphaBeta4[3]:$idAlphaBeta4[2]"})
        {
            $cellClonotype{$cellID}=$clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"};
            $clonotypeCount{$cellClonotype{$cellID}}++;
        }else
        {
            $cid++;
            my $nid=sprintf("C%04d",$cid);
            $cellClonotype{$cellID}=$nid;
            $clonotypeCount{$nid}=1;
            $clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"}=$nid;
            $clonotype{"$idAlphaBeta4[1]:$idAlphaBeta4[0]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"}=$nid;
            $clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[3]:$idAlphaBeta4[2]"}=$nid;
            $clonotype{"$idAlphaBeta4[1]:$idAlphaBeta4[0]:$idAlphaBeta4[3]:$idAlphaBeta4[2]"}=$nid;
        }
        push @inputData,$line;
    }
    foreach my $line (@inputData)
    {
        my @F=split /\t/,$line;
        my $cellID=$F[0];
        my $nid=$cellClonotype{$cellID};
        my $nid_count=$clonotypeCount{$nid};
        print $out join("\t", $line,"${opt_s}_$nid:$nid_count","$F[1]:$F[2]",$nid_count>1?"Clonal":"NoClonal",$F[2]>1?"Clonal":"NoClonal")."\n";
    }
}

sub usage
{
    die `pod2text $0`;
}
