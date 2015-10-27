#!/usr/bin/env perl
#============================================================================
# Name        		: report_gen.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Thu May 30 13:06:31 2013
# Last Modified By	: 
# Last Modified On	: Thu May 30 13:06:31 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl report_gen.pl [option] <outfile>
	-l	file list
	-t	template file
	-h	display this help and exit

=cut

use strict;
use warnings;
use utf8;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_l,$opt_t);
GetOptions("h"	=>\$opt_h,"l=s"=>\$opt_l,"t=s"=>\$opt_t);
if(@ARGV<1 || $opt_h) { usage(); }
if(!defined($opt_l)) { usage(); }
if(!defined($opt_t)) { usage(); }

### file list
my %h=();
open $in,$opt_l or die "$!";
while(<$in>)
{
	chomp;
	if(/^\s*$/) { next; }
	my ($id,$f)=split /=/;
	$h{$id}=$f;
}


my $outfile=shift @ARGV;
open $in,$opt_t or die "Cann't open file $opt_t ($!) \n";
open $out,">",$outfile or die "$!";
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^%%TABLE:(data_production)/)
	{
		#print "CAPTURED: $line\n";
		my $f=$h{$1};
		if($f)
		{
			my $_in;
			open $_in,$f or die "$!";
			my @d=<$_in>;
			chomp @d;
			@d=map { s/\t/ \& /g;$_ } @d;
			@d=map { s/\%/\\\%/g;$_ } @d;
			print $out "\\begin{tabular}{ccccccccc}\n";
			print $out "\\hline\n";
			print $out "$d[0] \\tabularnewline\n";
			print $out "\\hline\n";
			for(my $i=1;$i<=$#d;$i++) { print $out "$d[$i] \\tabularnewline\n"; }
			print $out "\\hline\n";
			print $out "\\end{tabular}\n";
			close $_in;
		}
	}elsif(/^%%TABLE:(mapping_stat)/)
	{
		my $f=$h{$1};
		if($f)
		{
			my $_in;
			open $_in,$f or die "$!";
			my @d=<$_in>;
			chomp @d;
			@d=map { s/\t/ \& /g;$_ } @d;
			@d=map { s/\%/\\\%/g;$_ } @d;
			print $out "\\begin{tabular}{lrrrr}\n";
			print $out "\\hline\n";
			print $out "$d[0] \\tabularnewline\n";
			print $out "\\hline\n";
			for(my $i=1;$i<=$#d;$i++) { print $out "$d[$i] \\tabularnewline\n"; }
			print $out "\\hline\n";
			print $out "\\end{tabular}\n";
			close $_in;
		}
	}
	elsif(/^%%TABLE:(snp_stat)/)
	{
		my $f=$h{$1};
		if($f)
		{
			my $_in;
			open $_in,$f or die "$!";
			my @d=<$_in>;
			chomp @d;
			@d=map { s/.+?\t//;$_ } @d;
			@d=map { s/\t/ \& /g;$_ } @d;
			@d=map { s/\%/\\\%/g;$_ } @d;
			print $out "\\begin{tabular}{c|c|c|cccc}\n";
			print $out "\\hline\n";
			print $out " & & & $d[0] \\tabularnewline\n";
			print $out "\\hline\n";
			print $out "\\multirow{28}{*}{Whole genome} & \\multirow{16}{*}{Region} & exonic & $d[1] \\tabularnewline\n";
			print $out " &  &  exonic,splicing & $d[2] \\tabularnewline\n";
			print $out " &  &  splicing & $d[3] \\tabularnewline\n";
			print $out " &  &  ncRNA exonic & $d[4] \\tabularnewline\n";
			print $out " &  &  ncRNA splicing & $d[5] \\tabularnewline\n";
			print $out " &  &  ncRNA UTR3 & $d[6] \\tabularnewline\n";
			print $out " &  &  ncRNA UTR5 & $d[7] \\tabularnewline\n";
			print $out " &  &  ncRNA intronic & $d[8] \\tabularnewline\n";
			print $out " &  &  UTR5 & $d[9] \\tabularnewline\n";
			print $out " &  &  UTR3 & $d[10] \\tabularnewline\n";
			print $out " &  &  UTR5,UTR3 & $d[11] \\tabularnewline\n";
			print $out " &  &  intronic & $d[12] \\tabularnewline\n";
			print $out " &  &  upstream & $d[13] \\tabularnewline\n";
			print $out " &  &  downstream & $d[14] \\tabularnewline\n";
			print $out " &  &  upstream,downstream & $d[15] \\tabularnewline\n";
			print $out " &  &  intergenic & $d[16] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & Total &  &  $d[17] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{3}{*}{Heterozygosity} &  Het & $d[18] \\tabularnewline\n";
			print $out " &  &  Hom & $d[19] \\tabularnewline\n";
			print $out " &  &  Het ratio & $d[20] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{3}{*}{Ts \\& Tv} &  transition & $d[21] \\tabularnewline\n";
			print $out " &  &  transvertion & $d[22] \\tabularnewline\n";
			print $out " &  &  ts/tv & $d[23] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{5}{*}{Novelty} &  in dbSNP (percentage) & $d[24] \\tabularnewline\n";
			print $out " &  &  novel & $d[25] \\tabularnewline\n";
			print $out " &  &  novel ts & $d[26] \\tabularnewline\n";
			print $out " &  &  novel tv & $d[27] \\tabularnewline\n";
			print $out " &  &  novel ts/tv & $d[28] \\tabularnewline\n";
			print $out "\\hline\n";
			print $out "\\multirow{17}{*}{Only CDS} & \\multirow{5}{*}{Function} & stopgain SNV &  $d[29] \\tabularnewline\n";
			print $out " &  &  stoploss SNV & $d[30] \\tabularnewline\n";
			print $out " &  &  nonsynonymous SNV & $d[31] \\tabularnewline\n";
			print $out " &  &  synonymous SNV & $d[32] \\tabularnewline\n";
			print $out " &  &  unknown & $d[33] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & Total &  &  $d[34] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{3}{*}{Heterozygosity} &  Het & $d[35] \\tabularnewline\n";
			print $out " &  &  Hom & $d[36] \\tabularnewline\n";
			print $out " &  &  Het ratio & $d[37] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{3}{*}{Ts \\& Tv} &  transition & $d[38] \\tabularnewline\n";
			print $out " &  &  transvertion & $d[39] \\tabularnewline\n";
			print $out " &  &  ts/tv & $d[40] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{5}{*}{Novelty} &  in dbSNP (percentage) & $d[41] \\tabularnewline\n";
			print $out " &  &  novel & $d[42] \\tabularnewline\n";
			print $out " &  &  novel ts & $d[43] \\tabularnewline\n";
			print $out " &  &  novel tv & $d[44] \\tabularnewline\n";
			print $out "&  &  novel ts/tv & $d[45] \\tabularnewline\n";
			print $out "\\hline\n";
			print $out "\\end{tabular}\n";
		}

	}elsif(/^%%TABLE:(indel_stat)/)
	{
		my $f=$h{$1};
		if($f)
		{
			my $_in;
			open $_in,$f or die "$!";
			my @d=<$_in>;
			chomp @d;
			@d=map { s/.+?\t//;$_ } @d;
			@d=map { s/\t/ \& /g;$_ } @d;
			@d=map { s/\%/\\\%/g;$_ } @d;
			print $out "\\begin{tabular}{c|c|c|cccc}\n";
			print $out "\\hline\n";
			print $out " & & & $d[0] \\tabularnewline\n";
			print $out "\\hline\n";
			print $out "\\multirow{28}{*}{Whole genome} & \\multirow{16}{*}{Region} & exonic & $d[1] \\tabularnewline\n";
			print $out " &  &  exonic,splicing & $d[2] \\tabularnewline\n";
			print $out " &  &  splicing & $d[3] \\tabularnewline\n";
			print $out " &  &  ncRNA exonic & $d[4] \\tabularnewline\n";
			print $out " &  &  ncRNA splicing & $d[5] \\tabularnewline\n";
			print $out " &  &  ncRNA UTR3 & $d[6] \\tabularnewline\n";
			print $out " &  &  ncRNA UTR5 & $d[7] \\tabularnewline\n";
			print $out " &  &  ncRNA intronic & $d[8] \\tabularnewline\n";
			print $out " &  &  UTR5 & $d[9] \\tabularnewline\n";
			print $out " &  &  UTR3 & $d[10] \\tabularnewline\n";
			print $out " &  &  UTR5,UTR3 & $d[11] \\tabularnewline\n";
			print $out " &  &  intronic & $d[12] \\tabularnewline\n";
			print $out " &  &  upstream & $d[13] \\tabularnewline\n";
			print $out " &  &  downstream & $d[14] \\tabularnewline\n";
			print $out " &  &  upstream,downstream & $d[15] \\tabularnewline\n";
			print $out " &  &  intergenic & $d[16] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & Total &  &  $d[17] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{3}{*}{Heterozygosity} & Het & $d[18] \\tabularnewline\n";
			print $out " & & Hom & $d[19] \\tabularnewline\n";
			print $out " & & Het ratio & $d[20] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{2}{*}{Novelty} &  in dbSNP (percentage) & $d[21] \\tabularnewline\n";
			print $out " & & novel & $d[22] \\tabularnewline\n";
			print $out "\\hline\n";
			print $out "\\multirow{17}{*}{Only CDS} & \\multirow{5}{*}{Function} & frameshift insertion & $d[23] \\tabularnewline\n";
			print $out " & & frameshift deletion & $d[24] \\tabularnewline\n";
			print $out " & & frameshift substitution & $d[25] \\tabularnewline\n";
			print $out " & & stopgain & $d[26] \\tabularnewline\n";
			print $out " & & stoploss & $d[27] \\tabularnewline\n";
			print $out " & & nonframeshift insertion & $d[28] \\tabularnewline\n";
			print $out " & & nonframeshift deletion & $d[29] \\tabularnewline\n";
			print $out " & & nonframeshift substitution & $d[30]  \\tabularnewline\n";
			print $out " & & unknown & $d[31] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & Total & & $d[32] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{3}{*}{Heterozygosity} & Het & $d[33] \\tabularnewline\n";
			print $out " & & Hom & $d[34] \\tabularnewline\n";
			print $out " & & Het ratio & $d[35] \\tabularnewline\n";
			print $out "\\cline{2-7}\n";
			print $out " & \\multirow{2}{*}{Novelty} & in dbSNP (percentage) & $d[36] \\tabularnewline\n";
			print $out " & & novel & $d[37] \\tabularnewline\n";
			print $out "\\hline\n";
			print $out "\\end{tabular}\n";
		}
	}
	else
	{
		print $out "$line\n"
	}
}
###################################################################



sub readList
{
	my $in;
	my ($pList,$infile)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		my $line=$_;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
	}
}

sub usage
{
	die `pod2text $0`;
}
