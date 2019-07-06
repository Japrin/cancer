#!/usr/bin/env perl
#============================================================================
# Name        		: refilter.expGene.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Thu Jun 25 09:20:06 2015
# Last Modified By	: 
# Last Modified On	: Thu Jun 25 09:20:06 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl refilter.expGene.pl [option] <infile>

    -e  min exp threshold required [default 3]
    -c  cancer types required [default "luad,lusc"]
    -n  min number of donors with the mutation required [default 1]
    -f  min percentage of donors with the mutation required [default 0.00 ]
    -m  min immunogenicity score required [default "NA" (not filter by this)]
    -p  mutation not in this position; disable by "NA" [default "1,2,9" ]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_e,$opt_c,$opt_n,$opt_f,$opt_m,$opt_p);
GetOptions("h"	=>\$opt_h,"e=f"=>\$opt_e,"c=s"=>\$opt_c,"n=i"=>\$opt_n,"f=f"=>\$opt_f,"m=s"=>\$opt_m,"p=s"=>\$opt_p);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_e)) { $opt_e=3; }
if(!defined($opt_c)) { $opt_c="luad,lusc"; }
if(!defined($opt_n)) { $opt_n=1; }
if(!defined($opt_f)) { $opt_f=0.00; }
if(!defined($opt_m)) { $opt_m="NA"; }
if(!defined($opt_p)) { $opt_p="1,2,9"; }
my @opt_c=split /,/,$opt_c;
my @opt_p=split /,/,$opt_p;
my %hCancerType=();
foreach (@opt_c) { $hCancerType{"\U$_"}=1;  }
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
print "wiltype peptide\twiltype logscore\twiltype affinity\twiltype Bind Level\tmutant peptide\tmutant logscore\tmutant affinity\tmutant Bind Level\tmutation pos in peptide\tpeptide length\tHLA allele\tchr\tpos\tref base\talt base\tgene symbol\ttranscript\tprotein\taa change\tmutation class\trecurrence\tcancer type\tindividual ID\tproject\tnValidExpFilter\texp filter\texpression\tOtherInfo1\tOtherInfo2\n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
	### 7: mutant IC50
	### 3: wild IC50
	### 9: mut pos in 9-mer
	### 21:number of total donor with the mutation
	### 25:number of donors with expressed mutant gene
	### 30:immunogenicity score

	#fmpwtpyRt       0.602   74      WB      fmpwtpyCt       0.624   58      WB      8       9       HLA-A02:01      1       82408790        C       T       LPHN2   ENST00000370717 ENSP00000359752 R179C   Missense_Mutation       4       skcm(3)|coadread(1)     skcm(TCGA-EE-A2MJ-06A-11D-A197-08;TCGA-EE-A3J5-06A-11D-A20D-08;TCGA-GN-A268-06A-11D-A196-08)|coadread(TCGA-AG-A002-01)  SAMPLE  3       PASS    skcm(TCGA-EE-A2MJ-06=8.22037042810397;TCGA-EE-A3J5-06=1.71058407582175;TCGA-GN-A268-06=16.2750655470734)|coadread(TCGA-AG-A002-01=12.7457022958738)     FMPWTPYCT       9       0.21203
	
	### parse record
	my $fieldTotalSample=$F[21];
	my @fieldTotalSample=split /\|/,$fieldTotalSample;
	my %hTotalSample=();
	my $nTotalSample=0;
	foreach (@fieldTotalSample)
	{
		my ($cType,$ss)=/(.+?)\((\d+)\)/;
		#print STDERR "$cType\t$ss\n";
		$cType="\U$cType";
		if($hCancerType{$cType})
		{
			$hTotalSample{$cType}=$ss;
			$nTotalSample+=$ss;
		}
	}
	my $fieldExpSample=$F[26];
	my @fieldExpSample=split /\|/,$fieldExpSample;
	my %hExpSample=();
	my $nExpSample=0;
	foreach (@fieldExpSample)
	{
		my ($cType,$ss)=/(.+?)\((.+?)\)/;
		$cType="\U$cType";
		if($hCancerType{$cType})
		{
			my @p=split /;/,$ss;
			foreach my $_p (@p)
			{
				my ($sample,$exp)=$_p=~/(.+)=(.+)/;
				if(!exists($hExpSample{$cType})) { $hExpSample{$cType}=[]; }
				push @{$hExpSample{$cType}},$exp;
				if($exp ne "NA" && $exp>=$opt_e)
				{
					$nExpSample++;
				}
			}
		}
	}
	my ($wildIC50,$mutantIC50,$mutPos,$immunoScore)=@F[2,6,8,29];
	### filtering procedure
	if(!($nTotalSample>=$opt_n && $nExpSample/$nTotalSample>$opt_f)) { next; }
	if($mutantIC50>=500) { next; }
	if($wildIC50<$mutantIC50) { next; }
	if($opt_p ne "NA")
	{
		my $_flag=0;
		foreach (@opt_p)
		{
			if($mutPos==$_) { $_flag=1; }
		}
		if($_flag) { next; }
	}
	#if($mutPos==1 || $mutPos==2 || $mutPos==9) { next; }
	if($opt_m ne "NA" && defined($immunoScore) && $immunoScore < $opt_m) { next; }
	
	splice @F,27,2;
	$line=join("\t",@F);
	# ' && $21>=3 && ( ($25/$21>0.2 && $30>0.05) || ($30>0.2) ) '
	print "$line\t$nTotalSample\t$nExpSample\n";
	###print "$line\t$nTotalSample\t$nExpSample\t$fieldTotalSample\t$fieldExpSample\n";

}
###################################################################




sub usage
{
    die `pod2text $0`;
}
