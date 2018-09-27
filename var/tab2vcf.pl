#!/usr/bin/env perl
#============================================================================
# Name        		: tab2vcf.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Thu Jan  1 19:10:05 2015
# Last Modified By	: 
# Last Modified On	: Thu Jan  1 19:10:05 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl tab2vcf.pl [option] <infile>
    
    -id	sample id(s), "," seperated [default "SAMPLEID"]
    -f  method related header
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Pod::Usage;
use FindBin qw($Bin);

my ($in,$in_f,$out);
my ($opt_h,$opt_f);
my ($id); 
GetOptions(
	"id=s"=>\$id,
	"h"=>\$opt_h,
	"f=s"=>\$opt_f
);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($id)) { $id="SAMPLEID"; }
#if(!defined($opt_v)) { $opt_v="general"; }
#my $hFile="$Bin/vcfHeader/$opt_v.header";

my @id = split ",", $id;
my $id_length = @id;

if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
if(defined($opt_f))
{
	open $in_f,$opt_f or die "$!";
	while(<$in_f>) { print "$_"; }
}
OutputANNOVARHeader();

print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".join("\t",@id)."\n";
my (@headers,%headers);
while (<$in>){
	chomp;
	my @record = split "\t";
	if ($record[-1] eq "Shared") {
		pop @record;
	}
	if (/^CHROM/) {
		foreach (0 .. $#record) {
			#$record[$_] =~ s/\..+$//;
			$headers{$record[$_]} = $_;
			push @headers, $record[$_];
		}
		
	}
	else {
        #my $dbsnp = $record[$headers{"avsnp142"}];
		my $dbsnp = $record[$headers{"avsnp150"}];
		foreach (0..6) {
			if (($_ == 2 ) and ($dbsnp)) {
				print "$dbsnp\t";
				next;
			}
			print "$record[$_]\t";
		}
		my $infoend = $#headers - 1 - $id_length;
		my @information;
		foreach (7..$infoend){
			next if ($record[$_] eq ".");
			next if ($headers[$_] eq "INFO");
			$record[$_] =~ s/ /_/g;
			$record[$_] =~ s/;/,/g;
			if ($record[$_] =~ /=/){
				$record[$_] = "("."$record[$_]".")";
			}
			push @information, "$headers[$_]=$record[$_]";
		}
		my $information = join ";", @information;
		$information=~s/nonsynonymous/missense/g;
		print $information;
		my $info = $record[$headers{"INFO"}];
		if ($info) {
			print ";$info";
		}
		my $format = $record[$headers{"FORMAT"}];
		print "\t$format";
		foreach (@id) {
			print "\t$record[$headers{$_}]";
		}
		print "\n";
	}
}

###################################################################

sub OutputANNOVARHeader
{
	print <<HERE
##INFO=<ID=Func,Number=.,Type=String,Description="Genome functional element(exonic,splicing,intergenic,intronic,ncRNA,upstream,downstream,UTR3,UTR5)">
##INFO=<ID=Gene,Number=.,Type=String,Description="Gene name">
##INFO=<ID=GeneDetail,Number=.,Type=String,Description="GeneDetail">
##INFO=<ID=ExonicFunc,Number=1,Type=String,Description="exonic variant functoin(frameshift insertion, frameshift deletion, frameshift block substitution, stopgain, stoploss, nonframeshift insertion, nonframeshift deletion, nonframeshift block substitution, missense, synonymous, unknown)">
##INFO=<ID=AAChange,Number=.,Type=String,Description="Amino acid change">
##INFO=<ID=cytoband,Number=.,Type=String,Description="the approximate location of bands seen on Giemsa-stained chromosomes">
##INFO=<ID=genomicSuperDups,Number=.,Type=String,Description="Segmental duplications in genome">
##INFO=<ID=phastConsElements46way,Number=.,Type=String,Description="conserved elements produced by the phastCons program based on a whole-genome alignment of vertebrates">
##INFO=<ID=evofold,Number=.,Type=String,Description="Conserve RNA prediction">
##INFO=<ID=wgRna,Number=.,Type=String,Description="snoRNA and miRNA annotation">
##INFO=<ID=targetScanS,Number=.,Type=String,Description="TargetScan generated miRNA target site predictions">
##INFO=<ID=tfbsConsSites,Number=.,Type=String,Description="transcription factor binding sites conserved in the human/mouse/rat alignment, based on transfac Matrix Database (v7.0)">
##INFO=<ID=esp6500siv2_all,Number=1,Type=Float,Description="variant frequency in the ESP6500 population">
##INFO=<ID=1000g2015aug_all,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=1000g2015aug_afr,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=1000g2015aug_eas,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=1000g2015aug_eur,Number=1,Type=Float,Description="variant frequency in the 1000G population">
##INFO=<ID=avsnp150,Number=.,Type=String,Description="dbSNP version 150">
##INFO=<ID=SIFT_score,Number=.,Type=String,Description="SIFT score">
##INFO=<ID=SIFT_pred,Number=.,Type=String,Description="SIFT score">
##INFO=<ID=Polyphen2_HDIV_score,Number=.,Type=String,Description="Polyphen2_HDIV_score, for complex phenotypes">
##INFO=<ID=Polyphen2_HDIV_pred,Number=.,Type=String,Description="Polyphen2_HDIV_pred, for complex phenotypes">
##INFO=<ID=Polyphen2_HVAR_score,Number=.,Type=String,Description="Polyphen2_HVAR_score, for Mendelian phenotypes">
##INFO=<ID=Polyphen2_HVAR_pred,Number=.,Type=String,Description="Polyphen2_HVAR_pred, for Mendelian phenotypes">
##INFO=<ID=LRT_score,Number=.,Type=String,Description="LRT_score">
##INFO=<ID=LRT_pred,Number=.,Type=String,Description="LRT_pred">
##INFO=<ID=MutationTaster_score,Number=.,Type=String,Description="MutationTaster_score">
##INFO=<ID=MutationTaster_pred,Number=.,Type=String,Description="MutationTaster_pred">
##INFO=<ID=MutationAssessor_score,Number=.,Type=String,Description="MutationAssessor_score">
##INFO=<ID=MutationAssessor_pred,Number=.,Type=String,Description="MutationAssessor_pred">
##INFO=<ID=FATHMM_score,Number=.,Type=String,Description="FATHMM_score">
##INFO=<ID=FATHMM_pred,Number=.,Type=String,Description="FATHMM_pred">
##INFO=<ID=RadialSVM_score,Number=.,Type=String,Description="RadialSVM_score">
##INFO=<ID=RadialSVM_pred,Number=.,Type=String,Description="RadialSVM_pred">
##INFO=<ID=LR_score,Number=.,Type=String,Description="LR_score">
##INFO=<ID=LR_pred,Number=.,Type=String,Description="LR_pred">
##INFO=<ID=VEST3_score,Number=.,Type=String,Description="VEST3_score">
##INFO=<ID=CADD_raw,Number=.,Type=String,Description="CADD_raw">
##INFO=<ID=CADD_phred,Number=.,Type=String,Description="CADD_phred">
##INFO=<ID=GERP++_RS,Number=.,Type=String,Description="GERP++_RS">
##INFO=<ID=phyloP46way_placental,Number=.,Type=String,Description="phyloP46way_placental">
##INFO=<ID=phyloP100way_vertebrate,Number=.,Type=String,Description="phyloP100way_vertebrate">
##INFO=<ID=SiPhy_29way_logOdds,Number=.,Type=String,Description="SiPhy_29way_logOdds">
##INFO=<ID=gerp++elem,Number=.,Type=Float,Description="conserved genomic regions by GERP++">
##INFO=<ID=gerp++gt2,Number=.,Type=Float,Description="conserved genomic sites by GERP++">
##INFO=<ID=cosmic70,Number=.,Type=String,Description="description of cosmic">
HERE


}


sub usage
{
    die `pod2text $0`;
}
