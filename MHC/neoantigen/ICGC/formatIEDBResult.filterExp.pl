#!/usr/bin/env perl
#============================================================================
# Name        		: formatIEDBResult.filterExp.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Mon Jul  6 17:38:08 2015
# Last Modified By	: 
# Last Modified On	: Mon Jul  6 17:38:08 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl formatIEDBResult.filterExp.pl [option] <infile>

    -e  gene expression file [required]
    -g  genename to ensemble gid file [default: /DBS/DB_temp/zhangLab/ICGC/PANCANCER/data/gencode.v19.annotation.hs37d5_chr.Genename2ENSG ]
    -k  knownTID to ensemble gid file [default: /DBS/DB_temp/zhangLab/ICGC/PANCANCER/data/knownToEnsembl.txt ]
    -t  gene expression threshold,  above this flag as "PASS" [default: 1]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_g,$opt_k,$opt_e,$opt_t);
GetOptions("h"	=>\$opt_h,"g=s"=>\$opt_g,"k=s"=>\$opt_k,"e=s"=>\$opt_e,"t=f"=>\$opt_t);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_g)) { $opt_g="/DBS/DB_temp/zhangLab/ICGC/PANCANCER/data/gencode.v19.annotation.hs37d5_chr.Genename2ENSG"; }
if(!defined($opt_k)) { $opt_k="/DBS/DB_temp/zhangLab/ICGC/PANCANCER/data/knownToEnsembl.txt"; }
if(!defined($opt_e)) { usage(); }
if(!defined($opt_t)) { $opt_t=1; }

if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"gzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
my %ID_MAPPING=();
readGeneMapping(\%ID_MAPPING,$opt_g,$opt_k);
my %GENE_EXP=();
readGeneExpression(\%GENE_EXP,$opt_e);
my $h=<$in>;
chomp $h;
print "$h\tensemble_gene_id\tgene_expression\texpFilter\n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
	my ($geneSymbol,$knownTID)=@F[12,13];
	my $exp;
	my $ensGID=$ID_MAPPING{$geneSymbol};
	if(!defined($ensGID))
	{
		$ensGID=$ID_MAPPING{$knownTID};
	}
	if(defined($ensGID))
	{
		$exp=$GENE_EXP{$ensGID};
	}else
	{
		$ensGID="NA";
	}
	if(!defined($exp)) { $exp="NA"; }
	my $filter="FAILED";
	if($exp ne "NA" && $exp>=$opt_t) { $filter="PASS"; }
	print "$line\t$ensGID\t$exp\t$filter\n";
}
###################################################################

sub readGeneMapping
{
    my $in;
    my ($pList,$infile,$infile2)=@_;
    open $in,$infile or die "Cann't open file $infile ($!) \n";
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
		my ($gname,$ensGID)=@F[0,1];
		$pList->{$gname}=$ensGID;
    }
	### patch
	if($infile2)
	{
		open $in,$infile2 or die "Cann't open file $infile2 ($!) \n";
		while(<$in>)
		{
			chomp;
			my $line=$_;
			if(/^\s*$/ || /^#/) { next; }
			my @F=split /\t/;
			my ($gname,$ensGID)=@F[0,1];
			$pList->{$gname}=$ensGID;
		}

	}
}

sub readGeneExpression
{
    my $in;
    my ($pList,$infile)=@_;

	if(defined($infile))
	{
		if($infile=~/\.gz$/) { open $in,"gzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
		elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
		else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
	}else
	{
		open $in,"-" or die "$!";
	}
	my $h=<$in>;
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
		#feature 3e69ca1f-5be0-4355-81e4-88e8ea2e1e56
		#ENSG00000000419.8       486.618167008
		my ($ensGID,$exp)=@F[0,1];
		$ensGID=~s/\.\d+//;
		$pList->{$ensGID}=$exp;

    }
}

sub usage
{
    die `pod2text $0`;
}
