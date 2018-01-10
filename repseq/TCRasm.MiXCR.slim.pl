#!/usr/bin/env perl
#============================================================================
# Name        		: TCRasm.MiXCR.slim.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Sat Jun  3 09:47:03 2017
# Last Modified By	: 
# Last Modified On	: Sat Jun  3 09:47:03 2017
# Copyright   		: Copyright (C) 2017
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl TCRasm.MiXCR.slim.pl [option] <infile> <outfile>

    -s  sample [default SAMPLE]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Dumper;

my ($in,$out);
my ($opt_h,$opt_s);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s);
if(@ARGV<2 || $opt_h) { usage(); }
my $infile=shift @ARGV;
my $outfile=shift @ARGV;
if(!defined($opt_s)){ $opt_s="SAMPLE"; }

open $in,$infile or die "Cann't open file $infile ($!) \n";
open $out,">",$outfile or die "Cann't open file $outfile ($!) \n";
my %ret=();
        
my $shh=<$in>;
my $n_alpha=0;
my $n_beta=0;
my @VDJ_A=();
my @VDJ_B=();
my @CDR3_nt_A=();
my @CDR3_nt_B=();
my @CDR3_aa_A=();
my @CDR3_aa_B=();
my @clone_count_A=();
my @clone_count_B=();
my @clone_freq_A=();
my @clone_freq_B=();
while(<$in>){
    chomp;
    my $sline=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @sF=split /\t/;
    if($sF[8]=~/^TRAC/){
        $n_alpha++;
        if($n_alpha<=2){
            ## [$vdj,$cdr3_nt,$cdr3_aa]
            my $pL=processMiTCRLine(\@sF);
            push @VDJ_A,$pL->[0];
            push @CDR3_nt_A,$pL->[1];
            push @CDR3_aa_A,$pL->[2];
            push @clone_count_A,$pL->[3];
            push @clone_freq_A,$pL->[4];
         }
    }elsif($sF[8]=~/^TRBC/){
        $n_beta++;
        if($n_beta<=2){
            my $pL=processMiTCRLine(\@sF);
            push @VDJ_B,$pL->[0];
            push @CDR3_nt_B,$pL->[1];
            push @CDR3_aa_B,$pL->[2];
            push @clone_count_B,$pL->[3];
            push @clone_freq_B,$pL->[4];
        }
    }
}
# VDJ
my @VDJ=();
while(@VDJ_A<2){ push @VDJ_A,"NA"; }
while(@VDJ_B<2){ push @VDJ_B,"NA"; }
push @VDJ,@VDJ_A;
push @VDJ,@VDJ_B;
# CDR3_nt
my @CDR3_nt=();
while(@CDR3_nt_A<2){ push @CDR3_nt_A,"NA"; }
while(@CDR3_nt_B<2){ push @CDR3_nt_B,"NA"; }
push @CDR3_nt,@CDR3_nt_A;
push @CDR3_nt,@CDR3_nt_B;
# CDR3_aa
my @CDR3_aa=();
while(@CDR3_aa_A<2){ push @CDR3_aa_A,"NA"; }
while(@CDR3_aa_B<2){ push @CDR3_aa_B,"NA"; }
push @CDR3_aa,@CDR3_aa_A;
push @CDR3_aa,@CDR3_aa_B;
# clone_count
my @clone_count=();
while(@clone_count_A<2){ push @clone_count_A,"NA"; }
while(@clone_count_B<2){ push @clone_count_B,"NA"; }
push @clone_count,@clone_count_A;
push @clone_count,@clone_count_B;
# clone_freq
my @clone_freq=();
while(@clone_freq_A<2){ push @clone_freq_A,"NA"; }
while(@clone_freq_B<2){ push @clone_freq_B,"NA"; }
push @clone_freq,@clone_freq_A;
push @clone_freq,@clone_freq_B;

print $out "sample\tVDJ_A1\tVDJ_A2\tVDJ_B1\tVDJ_B2\tcount_A1\tcount_A2\tcount_B1\tcount_B2\tfreq_A1\tfreq_A2\tfreq_B1\tfreq_B2\tCDR3_aa_A1\tCDR3_aa_A2\tCDR3_aa_B1\tCDR3_aa_B2\tCDR3_nt_A1\tCDR3_nt_A2\tCDR3_nt_B1\tCDR3_nt_B2\n";
print $out join("\t",$opt_s,@VDJ,@clone_count,@clone_freq,@CDR3_aa,@CDR3_nt)."\n";


##############################


sub processMiTCRLine
{
    my ($pF)=@_;
    my $v=(split /,/,$pF->[5])[0];
    $v=~s/\(.+\)//;
    my $d=(split /,/,$pF->[6])[0];
    if(defined($d)){
        $d=~s/\(.+\)//;
    }
    if(!defined($d) || $d eq "") { $d="."; }
    my $j=(split /,/,$pF->[7])[0];
    $j=~s/\(.+\)//;
    my $vdj=join("|",$v,$d,$j);
    my ($cdr3_nt,$cdr3_aa,$clone_count,$clone_freq)=@{$pF}[23,32,1,2];
    
    my $ret=[$vdj,$cdr3_nt,$cdr3_aa,$clone_count,$clone_freq];
    return($ret);
}

sub usage
{
    die `pod2text $0`;
}
