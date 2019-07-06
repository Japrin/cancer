#!/usr/bin/env perl
#============================================================================
# Name        		: maf2individualANNOVAR.pl
# Author      		: zhengl liangtao(tao2013@gmail.com)
# Version     		: v1.00
# Created On  		: Wed Jun 17 14:35:09 2015
# Last Modified By	: 
# Last Modified On	: Wed Jun 17 14:35:09 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl maf2individualANNOVAR.pl [option] <infile>

    -o  outDir [default .]
    -i  index of "project,sample";0-based [default 58,59]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;

my ($in,$out);
my ($opt_h,$opt_o,$opt_i);
GetOptions("h"	=>\$opt_h,"o=s"=>\$opt_o,"i=s"=>\$opt_i);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_o)) { $opt_o="."; }
if(!defined($opt_i)) { $opt_i="58,59"; }
my @opt_i=split /,/,$opt_i;
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"gzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
my $header=<$in>;
my $preSample="";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
	my ($chr,$beg,$end,$vclass,$ref,$alt,$gene,$transcript,$cDNAchange,$aachange,$proj,$sampleID)=@F[2,3,4,6,8,10,0,15,19,21,@opt_i];
	### insertion
	if($ref eq "-" && $alt ne "-") { $end=$beg; }
	if($sampleID ne $preSample)
	{
		mkpath "$opt_o/$proj";
		open $out,">","$opt_o/$proj/$proj.$sampleID.annovar" or die "$!";
		$preSample=$sampleID;
		print "process ${proj}::${sampleID} now ...\n";
	}
	print $out "$chr\t$beg\t$end\t$ref\t$alt\t$gene\t$transcript\t$aachange\t$vclass\t1\t$proj(1)\t$proj($sampleID)\t$sampleID\n";
}
print "completed successfully!\n";
###################################################################




sub usage
{
    die `pod2text $0`;
}
