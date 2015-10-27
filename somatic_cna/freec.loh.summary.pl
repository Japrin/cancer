#!/usr/bin/env perl
#============================================================================
# Name        		: freec.cnv.summary.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Fri Feb 14 15:37:19 2014
# Last Modified By	: 
# Last Modified On	: Fri Feb 14 15:37:19 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl freec.cnv.summary.pl [option] <infile>

	-s	sample [default "SAMPLE"]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_a,$opt_s);
GetOptions("h"	=>\$opt_h,"a"=>\$opt_a,"s=s"=>\$opt_s);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
#my $infile=shift @ARGV;
if(!defined($opt_s)) { $opt_s="SAMPLE"; }

my %stat=();
#open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#1       16931554        17235428        0.333333        0.666667        1       16931376        17235504        6       gain    AAAABBNA      somatic NA
	my ($cn,$cn_type,$somatic_status)=@F[8,9,12];
	$stat{$somatic_status}->{"$cn"}->{'count'}++;
	#if($opt_a)
	#{
	#	$stat{$somatic_status}->{"$cn"}->{'length'}+=$F[7]-$F[6];
	#}else
	{
		$stat{$somatic_status}->{"$cn"}->{'length'}+=$F[2]-$F[1];
	}
}
print "sample\tsomatic status\tCN\tcount\tlength\n";
foreach my $status (keys %stat)
{
	foreach my $cn (sort { $a<=>$b } keys %{$stat{$status}})
	{
		my $count=$stat{$status}->{$cn}->{'count'};
		my $length=$stat{$status}->{$cn}->{'length'};
		print "$opt_s\t$status\t$cn\t$count\t$length\n";
	}
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
