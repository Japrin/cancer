#!/usr/bin/perl
#============================================================================
# Name        		: compare.somatic.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Tue Sep  6 22:44:34 2011
# Last Modified By	: 
# Last Modified On	: Tue Sep  6 22:44:34 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl compare.somatic.pl [option] <inVCF> <valList>

	--t	inVCF type, vcf or bed [default vcf]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_t);
GetOptions("h"	=>\$opt_h,"t=s"=>\$opt_t);
if(@ARGV<2 || $opt_h) { usage(); }
my $infile=shift @ARGV;
my $valfile=shift @ARGV;
if(!defined($opt_t)) { $opt_t="vcf"; }
if($opt_t ne "vcf" && $opt_t ne "bed"){ usage(); }

my %list=();
readList(\%list,$valfile);

my $NA=0;
my $right=0;
my $wrong=0;
my $dataErr=0;
my $total=0;
my $comparable=0;
my $valTotal=0;
my $FN=0;

open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
	chomp;
	if(/^\s*$/ || /^#/) { next; }
	#chr1    17135864        .       G       A       101.40  TruthSensitivityTranche99.00to99.90     AAChange=NM_014675:c.G1102A:p.E
	my @field=split /\t/;
	my ($chr,$pos,$ref,$alt);
	if($opt_t eq "vcf"){ ($chr,$pos,$ref,$alt)=@field[0,1,3,4]; }
	elsif($opt_t eq "bed") { ($chr,$pos,$ref,$alt)=@field[0,2,3,4]; }
	my $key="$chr:$pos";
	$total++;
	if(!exists($list{$key})) { $NA++; }
	else
	{
		$list{$key}->{'comparable'}=1;
		if($list{$key}->{'somatic'} !~ /Yes|No/i) { printf STDERR "unknown somatic status\t$chr\t$pos\t%s\n",$list{$key}->{'somatic'};$NA++;next; }
		if($list{$key}->{'ref'} ne $ref) { printf STDERR "Unmatch reference\t$chr\t$pos\tval:%s\tinVCF:%s\n",$list{$key}->{'ref'},$ref;$NA++;next; }
		if($list{$key}->{'alt'} ne $alt) { printf STDERR "Unmatch calling\t$chr\t$pos\tval:%s\tinVCF:%s\n",$list{$key}->{'alt'},$alt;$wrong++;next; }
		if($list{$key}->{'somatic'} =~ /Yes/i) { $right++; }
		elsif($list{$key}->{'somatic'} =~ /No/i) { $wrong++;printf STDERR "Wrong prediction:\t%s\n",$_;next; }
	}
}
foreach (keys %list)
{
	$valTotal++;
	if($list{$_}->{'comparable'}==1) { $comparable++; }
	if($list{$_}->{'comparable'}==0 && $list{$_}->{'somatic'} =~/Yes/i) { $FN++;printf STDERR "FN\t$_\tref:%s\talt:%s\n",$list{$_}->{'ref'},$list{$_}->{'alt'}; }
}
printf "TP\t$right\n";
printf "FP\t$wrong\n";
printf "Rate\t%4.2f\n",($right+$wrong>0)?$right/($right+$wrong):-1;
printf "NA\t$NA\n";
printf "#Prediction\t$total\n";
printf "------\n";
printf "ValTotal\t$valTotal\n";
printf "Comparable\t$comparable\n";
printf "FN\t$FN\n";



############################################################################
sub readList
{
	my ($in);
	my ($pList,$infile)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		# #Chr    Start   Ref     Obs     somatic
		#chr11   60461415        A       G       Yes
		my @field=split /\t/;
		my ($chr,$pos,$ref,$alt,$somatic)=@field[0,1,2,3,4];
		my $key="$chr:$pos";
		$pList->{$key}={'ref'=>$ref,'alt'=>$alt,'somatic'=>$somatic,'comparable'=>0};
	}
}
sub usage
{
	die `pod2text $0`;
}
