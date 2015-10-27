#!/usr/bin/perl
#============================================================================
# Name        		: MentoCaloGenMutation.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Sun Sep 11 10:16:34 2011
# Last Modified By	: 
# Last Modified On	: Sun Sep 11 10:16:34 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl MentoCaloGenMutation.pl [option] <siteDir> <randomLineFile>

	-h	display this help and exit

	Note: randomeLineFile must be sorted by "_file" column and "_line" column

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;
use Cwd qw(abs_path);

my ($in,$in_site,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<2 || $opt_h) { usage(); }
my $siteDir=shift @ARGV;
my $infile=shift @ARGV;
$siteDir=abs_path $siteDir;

my $siteFile="";
my $line=0;
my $siteLine="";
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
	chomp;
	if(/^\s*$/) { next; }
	#CinCpG.bed.gz   1015261 C       T       0
	my @field=split /\t/;
	my ($_file,$_l,$ref,$alt,$_i)=@field[0,1,2,3,4];
	if($siteFile ne "$_file")
	{
		printf STDERR "info\topen another file\torigin:$siteFile\tnow:$_file\n";
		$siteFile=$_file;
		$line=0;
		#open another file
		if($_file =~ /\.gz$/) { open $in_site,"bgzip -cd $siteDir/$_file | " or die "Cann't open file $_file ($!)\n"; }
		else { open $in_site,"$siteDir/$_file" or die "Cann't open file $_file ($!)\n"; }
	}
	while($line != $_l)
	{
		$siteLine=<$in_site>;
		if(!defined($siteLine)) { die "out of range:\tsiteFile:$siteFile\tline:$_l\n"; }
		$line++
	}
	#chr1    58953   58954
	my ($chr,$beg,$end)=split /\t/,$siteLine;
	my $pos=$beg+1;
	my $iLast=@field-1;
	printf "$chr\t$pos\t$pos\t$ref\t$alt\t$_i\t$_file\t$_l\t%s\n",join("\t",@field[5 .. $iLast]);
}



############################################################################
sub usage
{
	die `pod2text $0`;
}
