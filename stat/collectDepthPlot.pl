#!/usr/bin/env perl
#============================================================================
# Name        		: collectDepthPlot.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Tue May 28 20:32:11 2013
# Last Modified By	: 
# Last Modified On	: Tue May 28 20:32:11 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl collectDepthPlot.pl [option] <infile>

	-i	col of  pos [default 0]
	-j	col of depth [defult 1]
	-t	file with header [default not]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my ($in,$out);
my ($opt_h,$opt_i,$opt_j,$opt_t);
GetOptions("h"	=>\$opt_h,"i=i"=>\$opt_i,"j=i"=>\$opt_j,"t"=>\$opt_t);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_i)) { $opt_i=0; }
if(!defined($opt_j)) { $opt_j=1; }

my %h=();
my @m=();
my @IDList=();
while(<>)
{
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	push @IDList,$F[0];
	my $m=readList(\%h,$F[1],$F[0]);
	push @m,$m;
}
my $MM=max(@m);
printf "depth\t%s\n",join("\t",@IDList);
for(my $i=0;$i<=$MM;$i++)
{
	print "$i";
	foreach my $_ID (@IDList)
	{
		printf "\t%s",defined($h{$i}->{$_ID})?$h{$i}->{$_ID}:0;
	}
	print "\n";
}
###################################################################



sub readList
{
	my $in;
	my ($pList,$infile,$ID)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	my $m=0;
	if($opt_t){ <$in>; }
	while(<$in>)
	{
		chomp;
		my $line=$_;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
		$pList->{$F[$opt_i]}->{$ID}=$F[$opt_j];
		$m=$F[0];
	}
	return $m;
}

sub usage
{
	die `pod2text $0`;
}
