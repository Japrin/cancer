#!/usr/bin/env perl
#============================================================================
# Name        		: merge.multiple.maf.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Tue Mar 31 17:21:06 2015
# Last Modified By	: 
# Last Modified On	: Tue Mar 31 17:21:06 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl merge.multiple.maf.pl [option] <class1:infile1> <class2:infile2> <class3:infile3> ...
	-s	sample [default "SAMPLE"]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_i,$opt_s);
GetOptions("h"	=>\$opt_h,"i=s"=>\$opt_i,"s=s"=>\$opt_s);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }

my %M=();
foreach (@ARGV)
{
	my ($c,$ifile)=/^(.+?):(.+)$/;
	readList(\%M,$c,$ifile);
}

foreach my $k (keys %M)
{
	
	print $M{$k}->{'basic'}."\t";
	print $M{$k}->{'total'}."\t";
	my @ctype=keys %{$M{$k}->{'cancerType'}};
	my @o=();
	foreach my $c (@ctype)
	{
		push @o,$c."(".$M{$k}->{'cancerType'}->{$c}->{'total'}.")";
	}
	print join("|",@o)."\t";
	@o=();
	foreach my $c (@ctype)
	{
		push @o,$c."(".join(";",@{$M{$k}->{'cancerType'}->{$c}->{'sample'}}).")";
	}
	print join("|",@o)."\t";
	print "$opt_s\n";
}


###################################################################

sub readList
{
    my $in;
    my $pList=shift @_;
    my $c=shift @_;
    my $infile=shift @_;
    open $in,$infile or die "Cann't open file $infile ($!) \n";
    while(<$in>)
    {
		chomp;
		my $line=$_;
		if(/^\s*$/ || /^#/ || /^Hugo_Symbol/) { next; }
		my @F=split /\t/;
    	my ($chr,$beg,$end,$hugo_symbol,$enst_id,$aaChange,$varClass,$sampleID,$ref,$alt)=@F[4,5,6,0,39,47,8,15,35,36];
		if($ref eq "-" && $alt ne "-") { $end=$beg; };
		my $key="$chr:$beg:$end:$ref:$alt";
		$pList->{$key}->{'total'}++;
		push @{$pList->{$key}->{'sample'}},$sampleID;
		if(!exists($pList->{$key}->{'basic'})) { $pList->{$key}->{'basic'}=join("\t",$chr,$beg,$end,$ref,$alt,$hugo_symbol,$enst_id,$aaChange,$varClass); }
		$pList->{$key}->{'cancerType'}->{$c}->{'total'}++;
		push @{$pList->{$key}->{'cancerType'}->{$c}->{'sample'}},$sampleID;

    }
}

sub usage
{
    die `pod2text $0`;
}
