#!/usr/bin/env perl
#============================================================================
# Name        		: formatIEDBResult.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Mon Jun 29 10:51:30 2015
# Last Modified By	: 
# Last Modified On	: Mon Jun 29 10:51:30 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl formatIEDBResult.pl [option] <infile>

    -b  bed file [required]
    -s  sampleID [default "SAMPLE"]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_b,$opt_s,$opt_p);
GetOptions("h"	=>\$opt_h,"b=s"=>\$opt_b,"s=s"=>\$opt_s,"p=s"=>\$opt_p);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_b)) { usage(); }
if(!defined($opt_p)) { $opt_p="PROJ"; }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;


open $in,"awk -F\"\t\" -v OFS=\"\t\" 'NR>1{print \$2,\$3-1,\$4,\$0}' $infile | intersectBed -a stdin -b $opt_b -wo -r -f 1 |" or die "$!";
print "allele\tpeptide\tmethod\tpercentile_rank\tann_ic50\tann_rank\tsmm_ic50\tsmm_rank\tcomblib_sidney2008_score\tcomblib_sidney2008_rank\tnetmhcpan_ic50\tnetmhcpan_rank\tgene_symbol\ttranscript\tmut_pos_in_peptide\tpeptdie_length\tcDNAChange\taaChange\tchr\tbeg\tend\tref\talt\tproject\tsample\n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
	## 1       3528    3537    HLA-A*02:05     1       3529    3537    9       LLMAALYQV       netmhcpan       0.2     -       -       -       -       -       -       3.87    0.2     1       3528    3537    HTR1A.uc011cqt.3.Mutant0001.ccd4a24b-d8cc-4686-9dee-c98b0c5a8d21_Pep0008 2/9 cDNAChange: [C272T], aaChange: [P91L], varSite alleles: [T], mutType: [SINGLEBASE_SUBSTITUTION], genomeMut: [5:63257275-63257275:G:A]      9
	my @partIEDB=@F[3..18];
	splice @partIEDB,1,4;
	my $partBed=$F[22];
	my ($gid,$tid,$mutPosInPep,$pepLen,$cdnaChange,$aaChange,$chr,$beg,$end,$ref,$alt)=$partBed=~/^(.+?)\.(.+?)\.Mutant.+ (.+?)\/(.+?) cDNAChange: \[(.+?)\], aaChange: \[(.+?)\],.+genomeMut: \[(.+?):(.+?)-(.+?):(.+?):(.+?)\]/;
	print join("\t",@partIEDB,($gid,$tid,$mutPosInPep,$pepLen,$cdnaChange,$aaChange,$chr,$beg,$end,$ref,$alt),$opt_p,$opt_s)."\n";
}
###################################################################




sub usage
{
    die `pod2text $0`;
}
