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
    ## output format:
    ## chr1 CREST   DEL 79470164    79470331    .   .   .   
    if($eventType eq "del" || $eventType eq "tandem_dup" || $eventType=~/invers*/)
    {
        #del|mechanism|cluster id|number of supporting read pairs|number of supporting split reads|chr|deletion boundary 1|deletion boundary 2|deletion size|homology sizes|annotation of break points
        # gff format
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[7]\t.\t.\t.\tSVID=$opt_p.$id;SVType=$F[0];SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=na;Distance=na;HomologySize=$F[9];AnnBP=$F[10];Range1=na;Range2=na;TCHR=na;TSTART=na\n";
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[6]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=na;Distance=na;HomologySize=$F[9];AnnBP=$F[10];Range1=na;Range2=na;TCHR=$F[5];TSTART=$F[7]\n";
        print "$F[5]\tmeerkat\t$F[0]\t$F[7]\t$F[7]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=na;Distance=na;HomologySize=$F[9];AnnBP=$F[10];Range1=na;Range2=na;TCHR=$F[5];TSTART=$F[6]\n";
        # annovar format
        #print "$F[5]\t$F[6]\t$F[7]\t0\t0\t$opt_p.$id\t$line\n";
    }elsif($eventType eq "del_invers" || $eventType =~/ins[so][ud]/)
    {
        #del_inss/o*     mechanism   cluster id  number of supporting read pairs     number of supporting split reads    chr     range of deletion (2 col)   deletion size   chr (donor)     range of insertion (2 col)  insert size     distance of deletion and insertion  homology at break points    annotation of break points
        if($F[6]>$F[7]) { my $tmp=$F[6];$F[6]=$F[7];$F[7]=$tmp; }
        if($F[10]>$F[11]) { my $tmp=$F[10];$F[10]=$F[11];$F[11]=$tmp; }
        # gff format
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[7]\t.\t.\t.\tSVID=$opt_p.$id;SVType=$F[0];SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=$F[13];HomologySize=$F[14];AnnBP=$F[15];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=na;TSTART=na\n";
        print "$F[9]\tmeerkat\t$F[0]\t$F[10]\t$F[11]\t.\t.\t.\tSVID=$opt_p.$id;SVType=$F[0];SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=$F[13];HomologySize=$F[14];AnnBP=$F[15];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=na;TSTART=na\n";
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[6]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=$F[13];HomologySize=$F[14];AnnBP=$F[15];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=$F[5];TSTART=$F[7]\n";
        print "$F[5]\tmeerkat\t$F[0]\t$F[7]\t$F[7]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=$F[13];HomologySize=$F[14];AnnBP=$F[15];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=$F[5];TSTART=$F[6]\n";
        print "$F[9]\tmeerkat\t$F[0]\t$F[10]\t$F[10]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=$F[13];HomologySize=$F[14];AnnBP=$F[15];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=$F[9];TSTART=$F[11]\n";
        print "$F[9]\tmeerkat\t$F[0]\t$F[11]\t$F[11]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=$F[13];HomologySize=$F[14];AnnBP=$F[15];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=$F[9];TSTART=$F[10]\n";
        # annovar format
        #print "$F[5]\t$F[6]\t$F[7]\t0\t0\t$opt_p.$id\t$line\n";
        #print "$F[9]\t$F[10]\t$F[11]\t0\t0\t$opt_p.$id\t$line\n";
    }elsif($eventType =~/ins[so]$/)
    {
        if($F[6]>$F[7]) { my $tmp=$F[6];$F[6]=$F[7];$F[7]=$tmp; }
        if($F[10]>$F[11]) { my $tmp=$F[10];$F[10]=$F[11];$F[11]=$tmp; }
        #del_inss/o  mechanism   cluster id  number of supporting read pairs     number of supporting split reads    chr     range of deletion (2 col)   deletion size   chr of insertion donor  range of insertion (2 col)  insert size     homology at break points    annotation of break points
        # gff format
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[7]\t.\t.\t.\tSVID=$opt_p.$id;SVType=$F[0];SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=na;HomologySize=$F[13];AnnBP=$F[14];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=na;TSTART=na\n";
        print "$F[9]\tmeerkat\t$F[0]\t$F[10]\t$F[11]\t.\t.\t.\tSVID=$opt_p.$id;SVType=$F[0];SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=na;HomologySize=$F[13];AnnBP=$F[14];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=na;TSTART=na\n";
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[6]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=na;HomologySize=$F[13];AnnBP=$F[14];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=$F[5];TSTART=$F[7]\n";
        print "$F[5]\tmeerkat\t$F[0]\t$F[7]\t$F[7]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=na;HomologySize=$F[13];AnnBP=$F[14];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=$F[5];TSTART=$F[6]\n";
        print "$F[9]\tmeerkat\t$F[0]\t$F[10]\t$F[10]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=na;HomologySize=$F[13];AnnBP=$F[14];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=$F[9];TSTART=$F[11]\n";
        print "$F[9]\tmeerkat\t$F[0]\t$F[11]\t$F[11]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=$F[12];Distance=na;HomologySize=$F[13];AnnBP=$F[14];Range1=$F[5]:$F[6]-$F[7];Range2=$F[9]:$F[10]-$F[11];TCHR=$F[9];TSTART=$F[10]\n";
        # annovar format
        #print "$F[5]\t$F[6]\t$F[7]\t0\t0\t$opt_p.$id\t$line\n";
        #print "$F[9]\t$F[10]\t$F[11]\t0\t0\t$opt_p.$id\t$line\n";
    }elsif($eventType eq "transl_inter")
    {
        #transl_inter    mechanism   cluster id  number of supporting read pairs     number of supporting split reads    chr of 1st cluster  boundary of 1st cluster     orientation of 1st cluster  chr of 2nd cluster  boundary of 2nd cluster     orientation of 2nd cluster  homology at break points    annotation of break points
        # gff format
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[6]\t.\t.\t.\tSVID=$opt_p.$id;SVType=$F[0];SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=$F[7];Size1=na;Size2=na;Distance=na;HomologySize=$F[11];AnnBP=$F[12];Range1=na;Range2=na;TCHR=$F[8];TSTART=$F[9]\n";
        print "$F[8]\tmeerkat\t$F[0]\t$F[9]\t$F[9]\t.\t.\t.\tSVID=$opt_p.$id;SVType=$F[0];SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=$F[10];Size1=na;Size2=na;Distance=na;HomologySize=$F[11];AnnBP=$F[12];Range1=na;Range2=na;TCHR=$F[5];TSTART=$F[6]\n";
        # annovar format
        #print "$F[5]\t$F[6]\t$F[6]\t0\t0\t$opt_p.$id\t$line\n";
        #print "$F[8]\t$F[9]\t$F[9]\t0\t0\t$opt_p.$id\t$line\n";
    }elsif($eventType eq "del_ins")
    {
        #del_ins     mechanism   cluster id  number of supporting read pairs     number of supporting split reads    chr    range of deletion (2 col)    deletion size  -   -   -    homology sizes  annotation of break points
        # gff format
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[7]\t.\t.\t.\tSVID=$opt_p.$id;SVType=$F[0];SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=na;Distance=na;HomologySize=$F[12];AnnBP=$F[13];Range1=na;Range2=na;TCHR=na;TSTART=na\n";
        print "$F[5]\tmeerkat\t$F[0]\t$F[6]\t$F[6]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=na;Distance=na;HomologySize=$F[12];AnnBP=$F[13];Range1=na;Range2=na;TCHR=$F[5];TSTART=$F[7]\n";
        print "$F[5]\tmeerkat\t$F[0]\t$F[7]\t$F[7]\t.\t.\t.\tSVID=$opt_p.$id;SVType=breakpoint;SVMechanism=$F[1];N_RP=$F[3];N_SP=$F[4];Ori1=na;Size1=$F[8];Size2=na;Distance=na;HomologySize=$F[12];AnnBP=$F[13];Range1=na;Range2=na;TCHR=$F[5];TSTART=$F[6]\n";
        # annovar format
        #print "$F[5]\t$F[6]\t$F[7]\t0\t0\t$opt_p.$id\t$line\n";
    }
    else
    {
        print STDERR "ERROR\t$line\n";
    }
}
###################################################################




sub usage
{
    die `pod2text $0`;
}
