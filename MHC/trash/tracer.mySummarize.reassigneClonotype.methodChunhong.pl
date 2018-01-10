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

    perl tracer.mySummarize.reassigneClonotype.methodChunhong.pl [option] <infile> <outfile>

    -s  sampleID [default: SAMPLE]
    -r  reorder by "main alpha > secondary alpha > main beta > secondary beta" [default: OFF]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_s,$opt_r);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s,"r"=>\$opt_r);
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
    if(!$opt_r) { print $out "$header\tC_strict\tC_share\tC_strict_b\tC_share_b\tmainAlpha\tmainBeta\n"; }
    else { print $out "$header\tC_strict\tC_share\tC_strict_b\tC_share_b\n"; }
    my $cid=0;
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
        my $cellID=$F[0];
        my @idAlphaBeta4=@F[3,9,15,21];
        my @cdr3AlphaBeta4=@F[4,10,16,22];
        my @TPMAlphaBeta4=@F[5,11,17,23];
        my @productiveAlphaBeta4=@F[6,12,18,24];
        my @inFrameAlphaBeta4=@F[7,13,19,25];
        my @stopCodonAlphaBeta4=@F[8,14,20,26];
        my $mainAlpha=-1;
        my $mainBeta=-1;
        #### main alpha
        if($productiveAlphaBeta4[0] eq "True" && $productiveAlphaBeta4[1]=~/False|NA/){
            $mainAlpha=0;
        }elsif($productiveAlphaBeta4[0]=~/False|NA/ && $productiveAlphaBeta4[1] eq "True"){
            $mainAlpha=1;
        }elsif($productiveAlphaBeta4[0] eq "True" && $productiveAlphaBeta4[1] eq "True"){
            $mainAlpha=$TPMAlphaBeta4[0]>=$TPMAlphaBeta4[1]?0:1;
        }
        #### main beta
        if($productiveAlphaBeta4[2] eq "True" && $productiveAlphaBeta4[3]=~/False|NA/){
            $mainBeta=2;
        }elsif($productiveAlphaBeta4[2]=~/False|NA/ && $productiveAlphaBeta4[3] eq "True"){
            $mainBeta=3;
        }elsif($productiveAlphaBeta4[2] eq "True" && $productiveAlphaBeta4[3] eq "True"){
            $mainBeta=$TPMAlphaBeta4[2]>=$TPMAlphaBeta4[3]?2:3;
        }
        if($mainAlpha==-1 || $mainBeta==-1)
        {
            print STDERR "mainAlpha=$mainAlpha\tmainBeta=$mainBeta\t$line\n";
            next;
        }
        if($clonotype{"$idAlphaBeta4[$mainAlpha]:$idAlphaBeta4[$mainBeta]"})
        {
            $cellClonotype{$cellID}={'nid'=>$clonotype{"$idAlphaBeta4[$mainAlpha]:$idAlphaBeta4[$mainBeta]"},"mainAlpha"=>"Alpha".($mainAlpha+1),"mainBeta"=>"Beta".($mainBeta-1)};
            $clonotypeCount{$cellClonotype{$cellID}->{'nid'}}++;
        }else
        {
            $cid++;
            my $nid=sprintf("C%04d",$cid);
            $cellClonotype{$cellID}={'nid'=>$nid,"mainAlpha"=>"Alpha".($mainAlpha+1),"mainBeta"=>"Beta".($mainBeta-1)};
            $clonotypeCount{$nid}=1;
            $clonotype{"$idAlphaBeta4[$mainAlpha]:$idAlphaBeta4[$mainBeta]"}=$nid;
        }
        if($opt_r){
            my $iMainAlpha=$mainAlpha;
            my $iSecondaryAlpha=1-$iMainAlpha;
            my $iMainBeta=$mainBeta;
            my $iSecondaryBeta=5-$iMainBeta;
            @idAlphaBeta4=@idAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
            @cdr3AlphaBeta4=@cdr3AlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
            @TPMAlphaBeta4=@TPMAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
            @productiveAlphaBeta4=@productiveAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
            @inFrameAlphaBeta4=@inFrameAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
            @stopCodonAlphaBeta4=@stopCodonAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
            my @tt=();
            for(my $ii=0;$ii<4;$ii++){
                push @tt,($idAlphaBeta4[$ii],$cdr3AlphaBeta4[$ii],$TPMAlphaBeta4[$ii],$productiveAlphaBeta4[$ii],$inFrameAlphaBeta4[$ii],$stopCodonAlphaBeta4[$ii]);
            }
            @F[3..26]=@tt;
            $line=join("\t",@F);
        }
        push @inputData,$line;
    }
    foreach my $line (@inputData)
    {
        my @F=split /\t/,$line;
        my $cellID=$F[0];
        my $nid=$cellClonotype{$cellID}->{'nid'};
        my $strMainAlpha=$cellClonotype{$cellID}->{'mainAlpha'};
        my $strMainBeta=$cellClonotype{$cellID}->{'mainBeta'};
        my $nid_count=$clonotypeCount{$nid};
        if(!$opt_r) { print $out join("\t", $line,"${opt_s}_$nid:$nid_count","$F[1]:$F[2]",$nid_count>1?"Clonal":"NoClonal",$F[2]>1?"Clonal":"NoClonal",$strMainAlpha,$strMainBeta)."\n"; }
        else{ print $out join("\t", $line,"${opt_s}_$nid:$nid_count","$F[1]:$F[2]",$nid_count>1?"Clonal":"NoClonal",$F[2]>1?"Clonal":"NoClonal")."\n"; }
    }
}

sub usage
{
    die `pod2text $0`;
}
