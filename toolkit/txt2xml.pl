#!/usr/bin/env perl
#============================================================================
# Name        		: txt2xml.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Mar 16 00:28:58 2015
# Last Modified By	: 
# Last Modified On	: Mon Mar 16 00:28:58 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl txt2xml.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

use XML::LibXML;

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
my $header=<$in>;
chomp $header;
my @header=split /\t/,$header;
my $doc = XML::LibXML::Document->new('1.0', 'utf-8');
my $root = $doc->createElement("QstatSummary");

while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
	my $ajob=$doc->createElement("Job");
	$root->appendChild($ajob);
    my @F=split /\t/;
	for(my $i=0;$i<@F;$i++)
	{
		my $n=$doc->createElement($header[$i]);
		$n->appendTextNode(($F[$i] eq "NA")?"":$F[$i]);
		$ajob->appendChild($n);
	}
}
$doc->setDocumentElement($root);
#$doc->toFH(*STDOUT);
print $doc->toString(4);
###################################################################




sub usage
{
    die `pod2text $0`;
}
