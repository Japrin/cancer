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

###-a  TPM threshold of alpha [default: 10]
###-b  TPM threshold of beta [default: 15]

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_s,$opt_r,$opt_a,$opt_b);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s,"r"=>\$opt_r,"a=i"=>\$opt_a,"b=i"=>\$opt_b);
if(@ARGV<2 || $opt_h) { usage(); }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
##if(!defined($opt_a)) { $opt_a=10; }
##if(!defined($opt_b)) { $opt_b=15; }

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
    #if(!$opt_r) { 
    #    ##print $out "$header\tC_strict\tC_share\tC_strict_b\tC_share_b\tmainAlpha\tmainBeta\n"; 
    #    print $out "$header\tFoundMajor\tC_strict\tC_share\tC_strict_b\tC_share_b\tmainAlpha\tmainBeta\n";
    #}
    #else 
    { 
        print $out "$header\tC_strict\tC_share\tC_strict_b\tC_share_b\n"; 
        #print $out "$header\tFoundMajor\tC_strict\tC_share\tC_strict_b\tC_share_b\n";
    }
    my $cid=0;
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
        my $cellID=$F[0];
        my @idAlphaBeta4=@F[3,11,19,27];
        my @cdr3AlphaBeta4=@F[4,12,20,28];
        my @cdr3ntAlphaBeta4=@F[5,13,21,29];
        my @VDJAlphaBeta4=@F[6,14,22,30];
        my @TPMAlphaBeta4=@F[7,15,23,31];
        my @productiveAlphaBeta4=@F[8,16,24,32];
        my @inFrameAlphaBeta4=@F[9,17,25,33];
        my @stopCodonAlphaBeta4=@F[10,18,26,34];
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
        my $foundMajor=1;
        if($mainAlpha==-1 || $mainBeta==-1)
        {
            print STDERR "No major: mainAlpha=$mainAlpha\tmainBeta=$mainBeta\t$line\n";
            $mainAlpha=0;
            $mainBeta=2;
            $foundMajor=0;
            #next;
        }
        #### if one productive alpha and one productive beta is the same, belong to the same clone
        my @ID_A_forClone=();
        my @ID_B_forClone=();
        if($productiveAlphaBeta4[0] eq "True") { push @ID_A_forClone,$idAlphaBeta4[0]; }
        if($productiveAlphaBeta4[1] eq "True") { push @ID_A_forClone,$idAlphaBeta4[1]; }
        if($productiveAlphaBeta4[2] eq "True") { push @ID_B_forClone,$idAlphaBeta4[2]; }
        if($productiveAlphaBeta4[3] eq "True") { push @ID_B_forClone,$idAlphaBeta4[3]; }
        my $beCounted=0;
        my $beFound=0;
        foreach my $_id_A (@ID_A_forClone){
            foreach my $_id_B (@ID_B_forClone){
                if($clonotype{"$_id_A:$_id_B"}){
                    ###$cellClonotype{$cellID}={'nid'=>$clonotype{"$_id_A:$_id_B"}};
                    if($beCounted==0){
                        ### move here
                        $cellClonotype{$cellID}={'nid'=>$clonotype{"$_id_A:$_id_B"}};
                        $clonotypeCount{$cellClonotype{$cellID}->{'nid'}}++;
                        $beCounted=1;
                    }
                    $beFound=1;
                }else{
                }
            }
        }
        if(!$beFound){
            $cid++;
            my $nid=sprintf("C%06d",$cid);
            if(@ID_A_forClone>0 && @ID_B_forClone>0){
                $cellClonotype{$cellID}={'nid'=>$nid};
                $clonotypeCount{$nid}=1;
                foreach my $_id_A (@ID_A_forClone){
                    foreach my $_id_B (@ID_B_forClone){
                        $clonotype{"$_id_A:$_id_B"}=$nid;
                    }
                }
            }else{
                $cellClonotype{$cellID}={'nid'=>"NoPPair"};
                #$clonotypeCount{$nid}=1;
            }
        }
        #### 

        #### if the four ids are the same, belong to the same clone
#        if($clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"})
#        #    || $clonotype{"$idAlphaBeta4[1]:$idAlphaBeta4[0]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"}
#        #    || $clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[3]:$idAlphaBeta4[2]"}
#        #    || $clonotype{"$idAlphaBeta4[1]:$idAlphaBeta4[0]:$idAlphaBeta4[3]:$idAlphaBeta4[2]"})
#        {
#            $cellClonotype{$cellID}={'nid'=>$clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"},
#                                     "mainAlpha"=>"Alpha".($mainAlpha+1),"mainBeta"=>"Beta".($mainBeta-1)};
#            $clonotypeCount{$cellClonotype{$cellID}->{'nid'}}++;
#        }else
#        {
#            $cid++;
#            my $nid=sprintf("C%06d",$cid);
#            $cellClonotype{$cellID}={'nid'=>$nid,"mainAlpha"=>"Alpha".($mainAlpha+1),"mainBeta"=>"Beta".($mainBeta-1)};
#            $clonotypeCount{$nid}=1;
#            $clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"}=$nid;
#            $clonotype{"$idAlphaBeta4[1]:$idAlphaBeta4[0]:$idAlphaBeta4[2]:$idAlphaBeta4[3]"}=$nid;
#            $clonotype{"$idAlphaBeta4[0]:$idAlphaBeta4[1]:$idAlphaBeta4[3]:$idAlphaBeta4[2]"}=$nid;
#            $clonotype{"$idAlphaBeta4[1]:$idAlphaBeta4[0]:$idAlphaBeta4[3]:$idAlphaBeta4[2]"}=$nid;
#        }
#        
#        if($opt_r)
#        {
#            my ($iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta)=(0,1,2,3);
#            if($opt_r){
#                $iMainAlpha=$mainAlpha;
#                $iSecondaryAlpha=1-$iMainAlpha;
#                $iMainBeta=$mainBeta;
#                $iSecondaryBeta=5-$iMainBeta;
#            }
#            @idAlphaBeta4=@idAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
#            @cdr3AlphaBeta4=@cdr3AlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
#            @cdr3ntAlphaBeta4=@cdr3ntAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
#            @VDJAlphaBeta4=@VDJAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
#            @TPMAlphaBeta4=@TPMAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
#            @productiveAlphaBeta4=@productiveAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
#            @inFrameAlphaBeta4=@inFrameAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
#            @stopCodonAlphaBeta4=@stopCodonAlphaBeta4[$iMainAlpha,$iSecondaryAlpha,$iMainBeta,$iSecondaryBeta];
#            my @tt=();
#            for(my $ii=0;$ii<4;$ii++){
#                push @tt,($idAlphaBeta4[$ii],$cdr3AlphaBeta4[$ii],$cdr3ntAlphaBeta4[$ii],$VDJAlphaBeta4[$ii],
#                    $TPMAlphaBeta4[$ii],$productiveAlphaBeta4[$ii],$inFrameAlphaBeta4[$ii],$stopCodonAlphaBeta4[$ii]);
#            }
#            @F[3..34]=@tt;
#            $line=join("\t",@F,$foundMajor?"TRUE":"FALSE");
#        }
        push @inputData,$line;
    }
    foreach my $line (@inputData)
    {
        my @F=split /\t/,$line;
        my $cellID=$F[0];
        my $nid="NA";
        my $strMainAlpha="NA";
        my $strMainBeta="NA";
        my $nid_count=0;
        $nid=$cellClonotype{$cellID}->{'nid'};
        #$strMainAlpha=$cellClonotype{$cellID}->{'mainAlpha'};
        #$strMainBeta=$cellClonotype{$cellID}->{'mainBeta'};
        $nid_count=$clonotypeCount{$nid};
        if(!defined($nid_count)){ $nid_count=0; }
        #if(!$opt_r) { 
        #    print $out join("\t", $line,"${opt_s}_$nid:$nid_count","$F[1]:$F[2]",$nid_count>1?"Clonal":"NoClonal",$F[2]>1?"Clonal":"NoClonal",$strMainAlpha,$strMainBeta)."\n"; 
        #    #print $out join("\t", $line)."\n"; 
        #}
        #else
        { 
            print $out join("\t", $line,"${opt_s}_$nid:$nid_count","$F[1]:$F[2]",$nid_count>1?"Clonal":"NoClonal",$F[2]>1?"Clonal":"NoClonal")."\n"; 
            #print $out join("\t", $line)."\n"; 
        }
    }
}

sub usage
{
    die `pod2text $0`;
}
