#!/usr/bin/env perl
#============================================================================
# Name        		: vcf_sort.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Fri Mar  8 17:13:59 2013
# Last Modified By	: 
# Last Modified On	: Fri Mar  8 17:13:59 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl vcf_sort.pl [option] <infile>

	-c	order file [if not set, order as: chrM, chr1, chr2,...]
	-i	chr column [default 0]
	-d	only chr from chrM,chr1,chr2,...chrX,chrY [just function when -c is OFF]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_c,$opt_i,$opt_d);
GetOptions("h"	=>\$opt_h,"c=s"=>\$opt_c,"i=i"=>\$opt_i,"d"=>\$opt_d);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;
if(!defined($opt_i)) { $opt_i=0; }

my @DEFAULT_ORDER=("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT", "GL000207.1", "GL000226.1", "GL000229.1", "GL000231.1", "GL000210.1", "GL000239.1", "GL000235.1", "GL000201.1", "GL000247.1", "GL000245.1", "GL000197.1", "GL000203.1", "GL000246.1", "GL000249.1", "GL000196.1", "GL000248.1", "GL000244.1", "GL000238.1", "GL000202.1", "GL000234.1", "GL000232.1", "GL000206.1", "GL000240.1", "GL000236.1", "GL000241.1", "GL000243.1", "GL000242.1", "GL000230.1", "GL000237.1", "GL000233.1", "GL000204.1", "GL000198.1", "GL000208.1", "GL000191.1", "GL000227.1", "GL000228.1", "GL000214.1", "GL000221.1", "GL000209.1", "GL000218.1", "GL000220.1", "GL000213.1", "GL000211.1", "GL000199.1", "GL000217.1", "GL000216.1", "GL000215.1", "GL000205.1", "GL000219.1", "GL000224.1", "GL000223.1", "GL000195.1", "GL000212.1", "GL000222.1", "GL000200.1", "GL000193.1", "GL000194.1", "GL000225.1", "GL000192.1", "NC_007605", "hs37d5");

my @DEFAULT_ORDER_D=("M", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y");

my @order=();
if(defined($opt_c))
{
	readList(\@order,$opt_c);
}else
{
	if($opt_d) { @order=@DEFAULT_ORDER_D; }
	else { @order=@DEFAULT_ORDER; }
}


my %h=();
my $TMPDir=dirname($infile);
my $r=int(rand(1e6));
my %tmpfile=();
open $in,$infile or die "$!";
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { print "$line\n"; next; }
	my @F=split /\t/;
	my $chr=$F[$opt_i];
	$chr =~ s/"//g;
	if(!exists($h{$chr})) 
	{ 
		my $tfile="$TMPDir/.tmp.$chr.$r";
		open $h{$chr},">",$tfile or die "$!";
		$tmpfile{$chr}=$tfile;
	}
	$out=$h{$chr};
	print $out "$line\n";
	#if($chr eq "chr7_gl000195_random") { print "##### OOO ####\t$line\n"; }
}
close $in;
foreach (keys %h)
{
	close $h{$_};
}
foreach my $chr (@order)
{
	if($tmpfile{$chr})
	{
		my $tfile=$tmpfile{$chr};
		open $in,$tfile or die "$!";
		while(<$in>)
		{
			print "$_";
		}
		close $in;
		##unlink $tfile;
	}
}
foreach (keys %tmpfile) { unlink $tmpfile{$_}; }


############################################################################
sub usage
{
	die `pod2text $0`;
}
sub readList
{
	my ($pList,$infile)=@_;
	my $in;
	open $in,$infile or die "$!";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
		push @$pList,$F[0];
	}
}


