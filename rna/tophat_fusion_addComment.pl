#!/usr/bin/perl
#============================================================================
# Name        		: tophat_fusion_addComment.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Nov 24 20:53:06 2014
# Last Modified By	: 
# Last Modified On	: Mon Nov 24 20:53:06 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl tophat_fusion_addComment.pl [option] <infile>
    -s	sample list
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

my %list=();
readList(\%list,$opt_s);

my %normalPanel=();

## first read, get normal panel
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    #120N    ANXA2   15      60639739        ANXA2P2 9       33625241        13      10      12      404.31
    my $sampleID=$F[0];
    my ($gene1,$gene2)=@F[1,4];
    my $cls=$list{$sampleID};
    if($cls eq "normal")
    {
    	$normalPanel{"$gene1-$gene2"}++;
	$normalPanel{"$gene2-$gene1"}++;
    }

}
## add comment
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    my $sampleID=$F[0];
    my ($gene1,$gene2)=@F[1,4];
    my $comment="";
    my $cls=$list{$sampleID};
    if($normalPanel{"$gene1-$gene2"} || $normalPanel{"$gene2-$gene1"})
    {
    	$comment="in normal panel";
    }else
    {
    	$comment="tumor specific";
    }
    print "$line\t$cls\t$comment\n";
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
	my ($key,$value)=@F[0,1];
	$pList->{$key}=$value;
    }
}

sub usage
{
    die `pod2text $0`;
}
