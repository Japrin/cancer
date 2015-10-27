#!/usr/bin/perl
#============================================================================
# Name        		: SplitSitesByContext.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Fri Sep  9 17:00:16 2011
# Last Modified By	: 
# Last Modified On	: Fri Sep  9 17:00:16 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl SplitSitesByContext.pl [option] <outDir> <infile(bed)> <fasta>
	
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;
use lib "/ifshk1/BC_CANCER/01bin/lib/perl/local/lib/perl5/site_perl";
use lib "/ifshk1/BC_CANCER/01bin/lib/perl/local/lib/perl5/x86_64-linux-thread-multi/"; 
use Bio::DB::Sam;
use Cwd qw(abs_path);

my ($in,$out,$out_A,$out_T,$out_C,$out_G,$out_CinCpG,$out_GinCpG,$out_CnotCpG,$out_GnotCpG);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<3 || $opt_h) { usage(); }
my $outDir=shift @ARGV;
my $infile=shift @ARGV;
my $fasta=shift @ARGV;

#A
#T
#C
#G
#CinCpG
#GinCpG
#CnotCpG
#GnotCpG
if(! -e $fasta) { usage(); }
mkpath $outDir;
$outDir=abs_path($outDir);
open $out_A,"| bgzip -c > $outDir/A.bed.gz" or die "Cann't open file ($!)\n"; 
open $out_T,"| bgzip -c > $outDir/T.bed.gz" or die "Cann't open file ($!)\n"; 
open $out_C,"| bgzip -c > $outDir/C.bed.gz" or die "Cann't open file ($!)\n"; 
open $out_G,"| bgzip -c > $outDir/G.bed.gz" or die "Cann't open file ($!)\n"; 
open $out_CinCpG,"| bgzip -c > $outDir/CinCpG.bed.gz" or die "Cann't open file ($!)\n"; 
open $out_GinCpG,"| bgzip -c > $outDir/GinCpG.bed.gz" or die "Cann't open file ($!)\n"; 
open $out_CnotCpG,"| bgzip -c > $outDir/CnotCpG.bed.gz" or die "Cann't open file ($!)\n"; 
open $out_GnotCpG,"| bgzip -c > $outDir/GnotCpG.bed.gz" or die "Cann't open file ($!)\n"; 
if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!)\n"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }

my $fai = Bio::DB::Sam::Fai->load($fasta);
while(<$in>)
{
	chomp;
	if(/^\s*$/) { next; }
	my @field=split /\t/;
	my ($chr,$beg,$end)=@field[0,1,2];
	my $pos=$beg+1;  ##bed format
	my $dna_string = $fai->fetch(sprintf("$chr:%d-%d",$pos-1,$pos+1));
	$dna_string="\U$dna_string";
	my ($preBase,$curBase,$nextBase)=split //,$dna_string;
	if($curBase eq "C")
	{
		if($nextBase eq "G") { printf $out_CinCpG "$_\n"; }
		else { printf $out_CnotCpG "$_\n"; }
	}elsif($curBase eq "G")
	{
		if($preBase eq "C") { printf $out_GinCpG "$_\n"; }
		else { printf $out_GnotCpG "$_\n"; }
	}
	if($curBase eq "A")
	{
		printf $out_A "$_\n"; 
	}elsif($curBase eq "T")
	{
		printf $out_T "$_\n"; 
	}elsif($curBase eq "C")
	{
		printf $out_C "$_\n"; 
	}elsif($curBase eq "G")
	{
		printf $out_G "$_\n"; 
	}

}

#`tabix -p bed -b 2 -e 3 $outDir/A.bed.gz`;
#`tabix -p bed -b 2 -e 3 $outDir/T.bed.gz`;
#`tabix -p bed -b 2 -e 3 $outDir/C.bed.gz`;
#`tabix -p bed -b 2 -e 3 $outDir/G.bed.gz`;
#`tabix -p bed -b 2 -e 3 $outDir/CinCpG.bed.gz`;
#`tabix -p bed -b 2 -e 3 $outDir/GinCpG.bed.gz`;
#`tabix -p bed -b 2 -e 3 $outDir/CnotCpG.bed.gz`;
#`tabix -p bed -b 2 -e 3 $outDir/GnotCpG.bed.gz`;


############################################################################
sub usage
{
	die `pod2text $0`;
}
