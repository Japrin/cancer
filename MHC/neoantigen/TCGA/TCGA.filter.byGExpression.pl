#!/usr/bin/env perl
#============================================================================
# Name        		: TCGA.filter.byGExpression.pl
# Author      		: zhengliangtao (tao2013@gmail.com)
# Version     		: v1.00
# Created On  		: Sat May  2 15:22:38 2015
# Last Modified By	: 
# Last Modified On	: Sat May  2 15:22:38 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl TCGA.filter.byGExpression.pl [option] <infile>
    -db  sqlite database of gene expression [required]
    -e   gene expression(TPM) threshold [default 3]
    -h   display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use DBI;

my ($in,$out);
my ($opt_h,$opt_db,$opt_e);
GetOptions("h"	=>\$opt_h,"db=s"=>\$opt_db,"e=f"=>\$opt_e);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_db)) { usage(); }
if(!defined($opt_e)) { $opt_e=3; }

if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}

my $driver   = "SQLite";
my $dsn = "DBI:$driver:dbname=$opt_db";
my $userid = "";
my $password = "";
my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 }) or die $DBI::errstr;
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
	my $geneSymbol=$F[15];
	#my @samples=$F[18]=~/(TCGA-.+?-.+?-.+?)-/g;
	#flsdSnmfi       0.882   3       SB      flsdYnmfi       0.915   2       SB      5       9       HLA-A02:01      7       48353926 	C       A       ABCA13  ENST00000435803 ENSP00000411096 S3260Y  Missense_Mutation       2       lusc(1)|kirc(1) lusc(TCGA-21-1070-01A-01D-1521-08)|kirc(TCGA-A3-3387-01A-01D-1534-10)   SAMPLE
	#
	my %TPM=();
	my @dat_cancerType=split /\|/,$F[22];
	my @cancerType=();
	my $filter="FAILED";
	my $nExpValidate=0;
	foreach (@dat_cancerType)
	{
		my ($cancerType,$dat_samples)=/(.+?)\((.+?)\)/;
		push @cancerType,$cancerType;
		my @samples=split /;/,$dat_samples;
		foreach (@samples)
		{
			my ($sampleID)=/(TCGA-.+?-.+?-..)/;
			my $stmt = qq(select tau from gene_rsem where sample == '$sampleID' AND gene_symbol == '$geneSymbol' limit 1;);
			my $sth = $dbh->prepare( $stmt );
			my $rv = $sth->execute() or die $DBI::errstr;
			if($rv < 0) { print $DBI::errstr; }
			my @row = $sth->fetchrow_array();
			if(defined($row[0])) 
			{ 
				my $_tpm=$row[0]*1e6;
				push @{$TPM{$cancerType}},"$sampleID=$_tpm"; 
				if($_tpm >= $opt_e) { $nExpValidate++; }
			}
			else { push @{$TPM{$cancerType}},"$sampleID=NA"; }
		}
	}
	if($nExpValidate) { $filter="PASS"; }
	my @o=();
	foreach my $cancer_type (@cancerType)
	{
		push @o,"$cancer_type(".join(";",@{$TPM{$cancer_type}}).")";
	}
	print join("\t",@F,$nExpValidate,$filter,join("|",@o))."\n";
}
$dbh->disconnect();
###################################################################




sub usage
{
    die `pod2text $0`;
}
