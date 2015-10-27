#!/usr/bin/perl
#============================================================================
# Name        		: breakdancer_chose_bed.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Mon Apr 15 16:09:16 2013
# Last Modified By	: 
# Last Modified On	: Mon Apr 15 16:09:16 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl breakdancer_chose_bed.pl [option] [infile]

	-ann	annotated sv file [required]
	-bed	bed file of alignemt [required]
	-gen	gene length file [default /WPS/GR/zhengliangtao/01bin/annovar/annovar_2013Feb25/humandb_hg19/hg19_refGene/hg19_refGene.gene.maxLength.txt]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_bed,$opt_ann,$opt_gen);
GetOptions("h"	=>\$opt_h,"bed=s"=>\$opt_bed,"ann=s"=>\$opt_ann,"gen=s"=>\$opt_gen);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_ann)) { usage(); }
if(!defined($opt_bed)) { usage(); }
if(!defined($opt_gen)) { $opt_gen="/WPS/GR/zhengliangtao/01bin/annovar/annovar_2013Feb25/humandb_hg19/hg19_refGene/hg19_refGene.gene.maxLength.txt"; }

my %list=();
readAnn(\%list,$opt_ann);
my %gen=();
readGen(\%gen,$opt_gen);

my @L=();
my @N=();
my @T=();
open $in,$opt_bed or die "$!";
while(1)
{
	my $line=<$in>;
	if(!$line) { last; }
	chomp $line;
	if($line=~/^track/)
	{
		if(@L>0)
		{
			### outut?
			if(exists($list{$N[-1]}))
			{
				my $g=$list{$N[-1]}->[0];
				my $size=$list{$N[-1]}->[1];
				my $svid=$list{$N[-1]}->[2];
				my $g_len=$gen{$g};
				if(!defined($g_len)) { $g_len=100000; printf STDERR "no gene length\t$g\n"; }
				$T[-1]="$T[-1]\tgene=$g\tgene_length=$g_len\tsv_size=$size\tSVID=$svid";
				print "$T[-1]\n";
				printf "%s\n",join("\n",@L);
			}
			#else { printf STDERR "no in ann\t$N[-1]\n"; }
			@L=();
		}
		#track name=chr1_869417_DEL_835  description="BreakDancer chr1 869417 DEL 835"   useScore=0
		my ($name)=$line=~/name=(.+?_.+?_.+?)_/;
		#print STDERR "$name\n";
		push @N,$name;
		push @T,$line;
	}else
	{
		push @L,$line;
	}
}
### output last one ?
if(@N && exists($list{$N[-1]}))
{
	my $g=$list{$N[-1]}->[0];
	my $size=$list{$N[-1]}->[1];
	my $svid=$list{$N[-1]}->[2];
	my $g_len=$gen{$g};
	if(!defined($g_len)) { $g_len=100000; printf STDERR "no gene length\t$g\n"; }
	$T[-1]="$T[-1]\tgene=$g\tgene_length=$g_len\tsv_size=$size\tSVID=$svid";
	print "$T[-1]\n";
	printf "%s\n",join("\n",@L);
}
#else { printf STDERR "no in ann\t$N[-1]\n"; }
@L=();


############################################################################
sub usage
{
	die `pod2text $0`;
}
sub readAnn
{
	my ($pList,$infile)=@_;
	my $in;
	open $in,$infile or die "$!";
	my $h=<$in>;
	while(<$in>)
	{
		chomp;
		my $line=$_;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
		my ($chr,$beg,$type,$gene,$region)=@F[0,1,6,3,4];
		$gene=~s/"//g;
		$gene=~s/\(.+//;
		$gene=(split /[;,]/,$gene)[0];
		#chr1    12898059        12898306        PRAMEF11(dist=6795),LOC649330(dist=8930)        intergenic      "23906:AluYa5(SINE)"    Deletion        .       Size 316;Support 6;SVID 1
		#
		my ($size) = /Size (.+?);/;
		my ($svid) = /SVID (\d+)/;
		if($type ne "Translocation")
		{
			if($region eq "intergenic" || $chr eq "chrM" || $chr eq "M") { next; }
			if($type eq "Deletion") { $type="DEL"; }
			elsif($type eq "Insertion") { $type="INS"; }
			elsif($type eq "Inversion") { $type="INV"; }
			else { printf STDERR "ERROR\t$line\n"; next; }
			my $name="${chr}_${beg}_${type}";
			#print "$name\t$line\n";
			$pList->{$name}=[$gene,$size,$svid];
		}else
		{
			my $line2=<$in>;
			chomp $line2;
			my @F2=split /\t/,$line2;
			my ($chr2,$beg2,$type2,$gene2,$region2)=@F2[0,1,6,3,4];
			if($region eq "intergenic" && $region2 eq "intergenic" || ($chr eq "chrM") || ($chr2 eq "M") ) { next; }
			if($chr eq $chr2){ $type="ITX"; }
			else { $type="CTX"; }
			my $name="${chr}_${beg}_${type}";
			my $name2="${chr2}_${beg2}_${type}";
			#if($beg eq "34354367") { printf STDERR "$name\t$name2\n$line\n$line2\n"; }
			#print "$name\t$line\n";
			#print "$name2\t$line2\n";
			$pList->{$name}=[$gene,$size,$svid];
			$pList->{$name2}=[$gene,$size,$svid];
		}
	}
}
sub readGen
{
	my ($pList,$infile)=@_;
	my $in;
	open $in,$infile or die "$!";
	while(<$in>)
	{
		chomp;
		my @F=split /\t/;
		$pList->{$F[1]}=$F[0];
	}
}
