#!/usr/bin/perl
#============================================================================
# Name        		: split.chrom.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Sat Apr 27 10:58:11 2013
# Last Modified By	: 
# Last Modified On	: Sat Apr 27 10:58:11 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl split.chrom.pl [option] [infile]

	-w	500000
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_w);
GetOptions("h"	=>\$opt_h,"w=i"=>\$opt_w);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_w)) { $opt_w=500000; }

while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	my ($chr,$len)=@F[0,1];
	for(my $i=0;$i<$len;)
	{
		my $j=$i+$opt_w;
		if($j>$len)
		{
			$j=$len;
			printf "$chr\t$i\t$j\n";
			last;
		}else
		{
			printf "$chr\t$i\t$j\n";
		}
		$i+=$opt_w;
	}
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
