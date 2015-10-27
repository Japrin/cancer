#!/usr/bin/perl
#============================================================================
# Name        		: vcf.info.inlist.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Mon Dec 16 20:17:23 2013
# Last Modified By	: 
# Last Modified On	: Mon Dec 16 20:17:23 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl vcf.info.inlist.pl [option] <infile> <file2>

	-i	info field [default: Gene]
	-j	index in file2 [default: 0]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_i,$opt_j);
GetOptions("h"	=>\$opt_h,"i=s"=>\$opt_i,"j=s"=>\$opt_j);
if(@ARGV<2 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
my $listFile2=shift @ARGV;

if(!defined($opt_i)) { $opt_i="Gene"; }
if(!defined($opt_j)) { $opt_j="0"; }

my %list2=();
my @arg_j=split(/,/,$opt_j);
readList(\%list2,$listFile2,@arg_j);

if($infile =~ /\.vcf\.gz$/) { open $in,"bgzip -cd $infile | " or die "$!"; }
else { open $in,$infile or die "$!"; }
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { print "$line\n"; next; }
	my @F=split /\t/;
	#1       8415555 .       C       G       .       PASS    Func=exonic;Gene=RERE;ExonicFunc=missense_SNV;AAChange=RERE:NM_001042682:exon12:c.G2929C:p.A977P,RERE:NM_001042681:exon22:c.G4591C:p.A1531P,RERE:NM_012102:exon23:c.G4591C:p.A1531P;evofold=(Score=69,Name=8693_0_+_69);phastConsElements46way=(Score=591,Name=lod=336);cytoband=1p36.23;avsift=0;ljb_pp2=0.791335;ljb2_pp2hdiv=1.0;ljb2_pp2hvar=0.999  GT:GQ:DP:RD:AD:FREQ:DP4 0/0:.:71:68:0:0%:35,33,0,0      0/1:.:79:43:34:44.16%:22,21,19,15
	if($F[7]=~/$opt_i=(.+?);/ || $F[7]=~/$opt_i=([^;]+)/)
	{	
		#print STDERR "$1\n";
		if($list2{$1})
		{
			print "$line\n";
		}
	}
	#else
	#{
	#	print STDERR "$line\n";
	#}
}
###################################################################

sub readList
{
	#my ($list,$infile,$index)=@_;
	my $list=shift @_;
	my $infile=shift @_;
	my @index=@_;
	my $in;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		my @field=split /\t/;
		my $malformed=0;
		for(my $i=0;$i<@index;$i++) { if(!defined($field[$index[$i]])) { $malformed=1; } }
		if($malformed) { next; }
		my $key=join(":",@field[@index]);
		#print "$key\n";
		#my $key=(split)[$index];
		$list->{$key}=$_;
	}
}


sub usage
{
	die `pod2text $0`;
}
