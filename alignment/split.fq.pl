#!/usr/bin/perl
#============================================================================
# Name        		: split.fq.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Mon May 20 11:23:25 2013
# Last Modified By	: 
# Last Modified On	: Mon May 20 11:23:25 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl split.fq.pl [option] <infile, fq or fq.gz> <output prefix>

	-n	number of reads per file [default 40000000]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_n);
GetOptions("h"	=>\$opt_h,"n=i"=>\$opt_n);
if(@ARGV<2 || $opt_h) { usage(); }
my $infile=shift @ARGV;
my $outPrefix=shift @ARGV;
if(!defined($opt_n)) { $opt_n=40000000; }

if($infile=~/\.fq\.gz$/) { open $in,"gzip -cd $infile | " or die "$!"; }
elsif($infile=~/\.fq$/) { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
else { usage(); }
my $c=0;
my $i=0;
while(<$in>)
{
	chomp;
	my $id=$_;
	my $seq=<$in>;
	chomp $seq;
	my $line3=<$in>;
	chomp $line3;
	my $qual=<$in>;
	chomp $qual;
	$c++;
	if($c % $opt_n ==1)
	{
		$i++;
		my $_ofile=sprintf "$outPrefix.S%06d.fq.gz",$i;
		open $out,"| gzip -c > $_ofile" or die "$!";
		print "$_ofile\n";
	}
	print $out "$id\n";
	print $out "$seq\n";
	print $out "$line3\n";
	print $out "$qual\n";
}
printf STDERR "Total reads:\t$c\n";
printf STDERR "Number reads per file:\t$opt_n\n";
printf STDERR "Number of files:\t$i\n";
###################################################################




sub usage
{
	die `pod2text $0`;
}
