#!/usr/bin/perl
#============================================================================
# Name        		: mbased.filterVPKE.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Fri Dec 19 16:05:20 2014
# Last Modified By	: 
# Last Modified On	: Fri Dec 19 16:05:20 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl mbased.filterVPKE.pl [option] <infile>
    -e	exon length file [required]
    -a	top percentile high VPKE gene will be filtered [default 0.02]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_e,$opt_a);
GetOptions("h"	=>\$opt_h,"e=s"=>\$opt_e,"a=f"=>\$opt_a);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_e)) { usage(); }
if(!defined($opt_a)) { $opt_a=0.02; }

my %e=();
readList(\%e,$opt_e);
my %m=();

## first time reading mut file
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    #1	898851	898852	C	T	KLHL17
    #1	906271	906272	A	C	PLEKHN1
    my ($chr,$beg,$end,$ref,$alt,$gene)=@F[0..5];
    $m{$gene}->{'mut'}++;
}
close $in;

my $totalGene=0;
foreach my $gene (keys %m)
{
	if(!defined($e{$gene}->{'length'}))
	{
		print STDERR "no length info available\t$gene\n";
		$m{$gene}->{'VPKE'}=$m{$gene}->{'mut'}*1000/1e9;
		next;
	}
	$m{$gene}->{'VPKE'}=$m{$gene}->{'mut'}*1000/$e{$gene}->{'length'};
	$totalGene++;
}
my $toFilter=int($totalGene*$opt_a);
my $i=0;
foreach my $gene (sort { $m{$b}->{'VPKE'} <=> $m{$a}->{'VPKE'} } keys %m)
{
	$i++;
	if($i<=$toFilter) { $m{$gene}->{'remove'}=1; }
	else { $m{$gene}->{'remove'}=0; }
}

## second time reading mut file
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    #1	898851	898852	C	T	KLHL17
    #1	906271	906272	A	C	PLEKHN1
    my ($chr,$beg,$end,$ref,$alt,$gene)=@F[0..5];
    my $VPKE=$m{$gene}->{'VPKE'};
    my $toRemove=$m{$gene}->{'remove'};
    if(!$toRemove)
    {
    	print join("\t",$chr,$beg,$end,$ref,$alt,$gene,$VPKE,$toRemove)."\n";
    }
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
	#A2ML1	12	6309
	my ($gene,$chr,$length)=@F[0..2];
	$pList->{$gene}->{'length'}=$length;
    }
}

sub usage
{
    die `pod2text $0`;
}
