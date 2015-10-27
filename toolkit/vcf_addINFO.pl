#!/usr/bin/perl
#============================================================================
# Name        		: vcf_addINFO.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Fri Mar  1 11:18:09 2013
# Last Modified By	: 
# Last Modified On	: Fri Mar  1 11:18:09 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl vcf_addINFO.pl [option] [in vcf file]

	-ID	info id [required]
	-l	info file [required]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_ID,$opt_l);
GetOptions("h"	=>\$opt_h,"ID=s"=>\$opt_ID,"l=s"=>\$opt_l);
if(@ARGV<0 || $opt_h || !defined($opt_ID) || !defined($opt_l)) { usage(); }

my %list=();
readList(\%list,$opt_l);

while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) 
	{
		if(/##INFO=<ID=DP,/)
		{
			printf "##INFO=<ID=$opt_ID,Number=.,Type=String,Description=\"\">\n";
		}
		print "$line\n";
		next; 
	}
	my @F=split /\t/;
	my $k="$F[0]:$F[1]";
	if(exists($list{$k}))
	{
		$F[7].=";$list{$k}";
	}
	printf "%s\n",join("\t",@F);
}



############################################################################
sub readList
{
	my ($pList,$infile)=@_;
	my $in;
	open $in,$infile or die "$!";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		#chr1    20945055        cosmic61=ID:COSM146436,OCCURENCE:1(stomach)
		my @F=split /\t/;
		my $key="$F[0]:$F[1]";
		$pList->{$key}=$F[2];
	}
}

sub usage
{
	die `pod2text $0`;
}
