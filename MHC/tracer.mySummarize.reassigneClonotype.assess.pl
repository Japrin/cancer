#!/usr/bin/env perl
#============================================================================
# Name        		: tracer.mySummarize.reassigneClonotype.assess.pl
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

    perl tracer.mySummarize.reassigneClonotype.assess.pl [option] <infile> <outfile>

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
    open $in,$infile or die "Cann't open file $infile ($!) \n";
    open $out,">","$outfile.fracAndRecoverRate" or die "$!";
    my $header=<$in>;
    chomp $header;
    my %Data=();
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
        my $sID=$F[0];
        my ($idA1,$idA2,$idB1,$idB2)=@F[3,9,15,21];
        my ($tpmA1,$tpmA2,$tpmB1,$tpmB2)=@F[5,11,17,23];
        if($sID!~/^(.+?)\.(rep.+?)\.D(.+)$/){
            print STDERR "cell name error: $line\n";
            next;
        }
        my ($cellID,$repID,$dFrac)=($1,$2,$3);
        ($idA1,$idA2)=sort {$b cmp $a} ($idA1,$idA2);
        ($idB1,$idB2)=sort {$b cmp $a} ($idB1,$idB2);
        ##print $out join("\t",$cellID,$cellID,$repID,$dFrac,$idA1,$idA2,$idB1,$idB2,$tpmA1,$tpmA2,$tpmB1,$tpmB2)."\n";
        $Data{$cellID}->{$dFrac}->{$repID}=[$idA1,$idA2,$idB1,$idB2,$tpmA1,$tpmA2,$tpmB1,$tpmB2];
    }
    close $in;

    print $out "cellID\tdFrac\tisConsistentA\tisConsistentB\tisRecoveredA\tisRecoveredB\n";
    foreach my $cellID (keys %Data)
    {
        my $mRecomA="";
        my $mRecomB="";
        foreach my $dFrac (sort { $b<=>$a } keys %{$Data{$cellID}})
        {
            my $recombA="";
            my $recombB="";
            my $isConsistentA=1;
            my $isConsistentB=1;
            my $isRecoveredA=0;
            my $isRecoveredB=0;
            foreach my $repID (sort keys %{$Data{$cellID}->{$dFrac}})
            {
                my $rA=join(":",@{$Data{$cellID}->{$dFrac}->{$repID}}[0,1]);
                my $rB=join(":",@{$Data{$cellID}->{$dFrac}->{$repID}}[2,3]);
                ###print "$cellID\t$dFrac\t$repID\t$rA\t$rB\n";
                if($recombA eq "") { $recombA=$rA; }
                if($recombB eq "") { $recombB=$rB; }
                if($rA ne $recombA){ $isConsistentA=0; }
                if($rB ne $recombB){ $isConsistentB=0; }
            }
            if($mRecomA eq "" && $isConsistentA){ $mRecomA=$recombA; }
            if($mRecomB eq "" && $isConsistentB){ $mRecomB=$recombB; }
            if($isConsistentA && $recombA eq $mRecomA){ $isRecoveredA=1; }
            if($isConsistentB && $recombB eq $mRecomB){ $isRecoveredB=1; }
            #print $out join("\t",$cellID,$dFrac,$isConsistentA?"TRUE":"FALSE",$isConsistentB?"TRUE":"FALSE",
            #    $isRecoveredA?"TRUE":"FALSE",$isRecoveredB?"TRUE":"FALSE")."\n";
            print $out join("\t",$cellID,$dFrac,$isConsistentA,$isConsistentB,$isRecoveredA,$isRecoveredB)."\n";
        }
    }
    close $out;

    open $out,">","$outfile.tpmAndMinDepth" or die "$!";
    print $out "CellID\tRecomb\tminFrac\tTPM\n";
    foreach my $cellID (keys %Data)
    {
        my @mRecom=("","","","");
        my @mTPM=("","","","");
        my @minFrac=(1,1,1,1);
        my @failed=(0,0,0,0);
        foreach my $dFrac (sort { $b<=>$a } keys %{$Data{$cellID}})
        {
            my @recomb=("","","","");
            my @isConsistent=(1,1,1,1);
            my @isRecovered=(0,0,0,0);
            my @tpm=(0,0,0,0);
            my $nRep=0;
            foreach my $repID (sort keys %{$Data{$cellID}->{$dFrac}})
            {
                for my $i (0..3)
                {
                    if($recomb[$i] eq "") { $recomb[$i]=$Data{$cellID}->{$dFrac}->{$repID}->[$i]; }
                    elsif($recomb[$i] ne $Data{$cellID}->{$dFrac}->{$repID}->[$i]){ $isConsistent[$i]=0; }
                    if("NA" ne $Data{$cellID}->{$dFrac}->{$repID}->[$i+4])
                    {
                        $tpm[$i]+=$Data{$cellID}->{$dFrac}->{$repID}->[$i+4];
                    }
                }
                $nRep+=1;
            }
            for my $i (0..3)
            {
                if($mRecom[$i] eq "" && $isConsistent[$i] && $recomb[$i] ne "NA")
                { 
                    $mRecom[$i]=$recomb[$i];
                    $mTPM[$i]=$tpm[$i]/$nRep; 
                }
                if($isConsistent[$i] && $recomb[$i] eq $mRecom[$i]){ $isRecovered[$i]=1; }
                if($isRecovered[$i] && $recomb[$i] ne "NA" && !$failed[$i]){
                    $minFrac[$i]=$dFrac;
                }else{
                    $failed[$i]=1;
                }
            }
        }
        for my $i (0..3){
            if($mRecom[$i] ne ""){
                print $out join("\t",$cellID,$mRecom[$i],$minFrac[$i],$mTPM[$i])."\n";
            }
        }
    }
    close $out;
    
}

sub usage
{
    die `pod2text $0`;
}
