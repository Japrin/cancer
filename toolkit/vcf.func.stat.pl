#!/usr/bin/env perl
#============================================================================
# Name        		: vcf.func.stat.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Fri Jul 13 16:31:59 2012
# Last Modified By	: 
# Last Modified On	: Fri Jul 13 16:31:59 2012
# Copyright   		: Copyright (C) 2012
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl vcf.func.stat.pl [option] <vcf file>

	-f	apply filter [default not]
	-s	sample column [default -1]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_f,$opt_s);
GetOptions("h"	=>\$opt_h,'f'=>\$opt_f,"s=i"=>\$opt_s);
if(@ARGV<1 || $opt_h) { usage(); }
if(!defined($opt_s)) { $opt_s=-1; }

my $vcfFile=shift @ARGV;

my %ttt=('AG'=>'ts',
	'GA'=>'ts',
	'TC'=>'ts',
	'CT'=>'ts',
	'AT'=>'tv',
	'TA'=>'tv',
	'AC'=>'tv',
	'CA'=>'tv',
	'GT'=>'tv',
	'TG'=>'tv',
	'GC'=>'tv',
	'CG'=>'tv',
);

my %func=();
my %exonFunc=();
my $funcTotal=0;
my $exonFuncTotal=0;
my $funcInDB=0;
my $exonFuncInDB=0;
my $funcHet=0;
my $funcHom=0;
my $exonFuncHet=0;
my $exonFuncHom=0;
my $funcTs=0;
my $funcTv=0;
my $exonFuncTs=0;
my $exonFuncTv=0;
my $novelFuncTs=0;
my $novelFuncTv=0;
my $novelExonFuncTs=0;
my $novelExonFuncTv=0;


my @funcCata=("exonic","exonic,splicing","splicing","ncRNA_exonic","ncRNA_splicing","ncRNA_UTR3","ncRNA_UTR5","ncRNA_intronic","UTR5","UTR3","UTR5,UTR3","intronic","upstream","downstream","upstream,downstream","intergenic");
#my @exonFuncCata=("frameshift_insertion","frameshift_deletion","frameshift_substitution","stopgain_SNV","stoploss_SNV","nonframeshift_insertion","nonframeshift_deletion","nonframeshift_substitution","missense_SNV","synonymous_SNV","unknown");
my @exonFuncCata=("frameshift_insertion","frameshift_deletion","frameshift_substitution","stopgain","stoploss","nonframeshift_insertion","nonframeshift_deletion","nonframeshift_substitution","missense_SNV","synonymous_SNV","unknown");

my $sampleID="";
if($vcfFile=~/\.gz$/) { open $in,"bgzip -cd $vcfFile | " or die "$!"; }
else { open $in,$vcfFile or die "$!"; }
while(<$in>)
{
	chomp;
	if(!/^#CHROM/) { next; }
	$sampleID=(split /\t/)[$opt_s];
	last;
}

while(<$in>)
{
	chomp;
	if(/^\s*$/) { next; }
	my @field=split /\t/;

	if($opt_f && $field[6] ne "." && $field[6] ne "PASS") { next; }
	
	my $inExon=0;
	if(/\bFunc=(.+?);/)
	{
		$func{$1}++;	
		$funcTotal++;
		if(/ExonicFunc=(.+?);/)
		{
			$inExon=1;
			my $_ef=$1;
			if($_ef eq "nonsynonymous_SNV") { $_ef="missense_SNV"; }
			$exonFunc{$_ef}++;
			$exonFuncTotal++;
		}
	}
	### het or hom
	if($field[$opt_s]=~/0\/1/ || /\shet\s/) 
	{ 
		$funcHet++;
		if($inExon) { $exonFuncHet++; }
	}
	else 
	{ 
		$funcHom++;
		if($inExon) { $exonFuncHom++; }
	}
	### ts or tv
	my $isTs=-1;
	my $ref=$field[3];
	my $alt=(split /,/,$field[4])[0];
	if(length($ref)==length($alt) && length($ref)==1 && $ref ne "-" && $alt ne "-")
	{
		my $cc=$ttt{"$ref$alt"};
		if($cc eq "ts") 
		{ 
			$isTs=1;
			$funcTs++; 
			if($inExon) { $exonFuncTs++; }
		}
		else 
		{ 
			$isTs=0;
			$funcTv++; 
			if($inExon) { $exonFuncTv++; }
		}
	}
	#else
	#{
	#	printf STDERR "$field[3]$field[4]\n";
	#}
	### in dbSNP or not
	if($field[2] ne ".")
	{
		$funcInDB++;
		if($inExon) { $exonFuncInDB++; }
	}else
	{
		if($isTs==1) 
		{ 
			$novelFuncTs++;
			if($inExon) { $novelExonFuncTs++; }
		}elsif($isTs==0)
		{
			$novelFuncTv++;
			if($inExon) { $novelExonFuncTv++; }
		}

	}
}

print "#genome\n";
print "#sampleID";
foreach my $f (@funcCata)
{
	print "\t$f";
}
print "\tTotal\tHet\tHom\tHet_ratio\ttransition\ttransvertion\tts/tv\tin dbSNP (percentage)\tnovel\tnovel ts\tnovel tv\tnovel ts/tv\n";

print "$sampleID";
foreach my $f (@funcCata)
{
	print "\t",$func{$f}?$func{$f}:0;
}
printf "\t$funcTotal\t$funcHet\t$funcHom\t%4.2f\t$funcTs\t$funcTv\t%4.2f\t$funcInDB (%4.2f%%)\t%d\t$novelFuncTs\t$novelFuncTv\t%4.2f\n",$funcHom==0?"-1":($funcHet/$funcTotal),$funcTv==0?"-1":($funcTs/$funcTv),$funcTotal==0?"-1":($funcInDB*100/$funcTotal),$funcTotal-$funcInDB,$novelFuncTv==0?"-1":($novelFuncTs/$novelFuncTv);

print "#exome\n";
print "#sampleID";
foreach my $f (@exonFuncCata)
{
	print "\t$f";
}
print "\tTotal\tHet\tHom\tHet_ratio\ttransition\ttransvertion\tts/tv\tin dbSNP (percentage)\tnovel\tnovel ts\tnovel tv\tnovel ts/tv\tNS:S ratio\n";

print "$sampleID";
foreach my $f (@exonFuncCata)
{
	if(!defined($exonFunc{$f})) { $exonFunc{$f}=0; }
	print "\t",$exonFunc{$f};
}
printf "\t$exonFuncTotal\t$exonFuncHet\t$exonFuncHom\t%4.2f\t$exonFuncTs\t$exonFuncTv\t%4.2f\t$exonFuncInDB (%4.2f%%)\t%d\t$novelExonFuncTs\t$novelExonFuncTv\t%4.2f\t%4.2f\n",$exonFuncHom==0?"-1":($exonFuncHet/$exonFuncTotal),$exonFuncTv==0?"-1":($exonFuncTs/$exonFuncTv),$exonFuncTotal==0?"-1":($exonFuncInDB*100/$exonFuncTotal),$exonFuncTotal-$exonFuncInDB,$novelExonFuncTv==0?"-1":($novelExonFuncTs/$novelExonFuncTv),$exonFunc{"synonymous_SNV"}==0?"-1":(($exonFunc{"missense_SNV"}+$exonFunc{"stopgain_SNV"}+$exonFunc{"stoploss_SNV"})/$exonFunc{"synonymous_SNV"});


############################################################################
sub usage
{
	die `pod2text $0`;
}
