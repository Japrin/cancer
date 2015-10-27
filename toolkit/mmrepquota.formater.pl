#!/usr/bin/env perl
#============================================================================
# Name        		: mmrepquota.formater.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Mar 16 09:04:55 2015
# Last Modified By	: 
# Last Modified On	: Mon Mar 16 09:04:55 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl mmrepquota.formater.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;


if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}

print "Device\tName\ttype\tBlock_current_usage-KB\tBlock_soft_limit-KB\tBlock_hard_limit-KB\tBlock_in_doubt\tBlock_grace_period\tFile_current_number\tFile_soft_limit\tFile_hard_limit\tFile_in_doubt\tFile_grace_period\tentryType\n";
$/="***";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
	if(/quotas on (.+)\b/)
	{
		process($line,$1);
	}
}
$/="\n";
###################################################################

sub process
{
	#my $_a=$/;
	my ($blockTxt,$dev)=@_;
	my @lines=split /\n/,$blockTxt;
	for my $line (@lines)
	{
		if($line!~/^Name/ && $line!~/Block Limits/ && $line!~/Report for/) 
		{
		   $line=~s/default on/default_on/;
		   $line=~s/default off/default_off/;
		   $line=~s/\|//g;
		   $line=~s/\s+/\t/g;
		   my @F=split /\t/,$line;
		   splice @F,0,0,$dev;
		   print join("\t",@F)."\n";
		}
	}	
	#$/=$_a;
}


sub usage
{
    die `pod2text $0`;
}
