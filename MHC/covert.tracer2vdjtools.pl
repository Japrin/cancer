#!/usr/bin/env perl
#============================================================================
# Name        		: covert2vdjtools.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Sat Feb 25 22:13:47 2017
# Last Modified By	: 
# Last Modified On	: Sat Feb 25 22:13:47 2017
# Copyright   		: Copyright (C) 2017
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl covert2vdjtools.pl [option] <infile>

    -b  only use beta chain
    -a  only use alpha chain
    -i  stratified by column [INT]; 0-based (default -1)
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_b,$opt_a,$opt_i);
GetOptions("h"	=>\$opt_h,"b"=>\$opt_b,"a"=>\$opt_a,"i=i"=>\$opt_i);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
my $gChainI=0;
if($opt_b) { $gChainI=2; }
my $gChainJ=4;
if($opt_a) { $gChainJ=2; }
if(!defined($opt_i)) { $opt_i=-1; }

if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
my $hh=<$in>;
my %D=();
my $nTotal=0;
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    my @ID4=@F[3,11,19,27];
    my @CDR3aa4=@F[4,12,20,28];
    my @CDR3nt4=@F[5,13,21,29];
    my @VDJ4=@F[6,14,22,30];
    my @prod4=@F[8,16,24,32];

    for(my $i=$gChainI; $i<$gChainJ; $i++){
        if($prod4[$i] eq "True" && $CDR3nt4[$i] ne "" && $CDR3aa4[$i]!~/Couldn/){
            if(exists($D{$ID4[$i]})){
                $D{$ID4[$i]}->{'count'}++;
            }else{
                my ($v_seg,$d_seg,$j_seg)=split /\|/,$VDJ4[$i];
                $v_seg=~s/\*.+//;
                $d_seg=~s/\*.+//;
                $j_seg=~s/\*.+//;
                if($d_seg eq "N/A"){ $d_seg="."; }
                $CDR3nt4[$i]=~s/.........$//;
                $CDR3aa4[$i]=~s/...$//;
                $D{$ID4[$i]}={'count'=>1,'cdr3nt'=>uc($CDR3nt4[$i]),'cdr3aa'=>$CDR3aa4[$i],
                    'v'=>$v_seg,'d'=>$d_seg,'j'=>$j_seg};
            }
            $nTotal++;
        }
    }
}
print "#count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\n";
for my $cid (sort { $D{$b}->{'count'} <=> $D{$a}->{'count'} } keys %D){
    print(join("\t",$D{$cid}->{'count'},$D{$cid}->{'count'}/$nTotal,$D{$cid}->{'cdr3nt'},$D{$cid}->{'cdr3aa'},$D{$cid}->{'v'},$D{$cid}->{'d'},$D{$cid}->{'j'})."\n");
}
###################################################################



sub readList
{
    my $in;
    my ($pList,$infile)=@_;
    open $in,$infile or die "Cann't open file $infile ($!) \n";
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
    }
}

sub usage
{
    die `pod2text $0`;
}
