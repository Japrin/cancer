#!/usr/bin/perl
#============================================================================
# Name        		: strelka.addFreq.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Wed Dec 31 20:17:12 2014
# Last Modified By	: 
# Last Modified On	: Wed Dec 31 20:17:12 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl strelka.addFreq.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
my $infile=shift @ARGV;

my $varType="";
if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/) { print "$line\n"; next; }
    if(/^#/)
    {
    	if(/^#CHROM\tPOS\tID/)
	{
		print "##INFO=<ID=NORMAL_FREQ,Number=1,Type=Float,Description=\"mutation frequency in normal\">\n";
		print "##INFO=<ID=TUMOR_FREQ,Number=1,Type=Float,Description=\"mutation frequency in tumor\">\n";
	}
    	print "$line\n"; 
	next;
    }
    my @F=split /\t/;
    my ($chr,$pos,$ref,$alt,$info,$sFormat,$v_normal,$v_tumor)=@F[0,1,3,4,7,8,9,10];
    my $freq_normal="NA";
    my $freq_tumor="NA";
    if($varType eq "SNV")
    {
	$freq_normal=getSNVFreq($ref,$alt,$v_normal);
	$freq_tumor =getSNVFreq($ref,$alt,$v_tumor);
	$info="$info;NORMAL_FREQ=$freq_normal;TUMOR_FREQ=$freq_tumor";
    }elsif($varType eq "INDEL")
    {
	$freq_normal=getINDELFreq($v_normal);
	$freq_tumor =getINDELFreq($v_tumor);
	$info="$info;NORMAL_FREQ=$freq_normal;TUMOR_FREQ=$freq_tumor";

    }else
    {
	if($sFormat=~/DP:FDP:SDP:SUBDP:AU:CU:GU:TU/) 
	{ 
		$varType="SNV"; 
		$freq_normal=getSNVFreq($ref,$alt,$v_normal);
		$freq_tumor =getSNVFreq($ref,$alt,$v_tumor);
		$info="$info;NORMAL_FREQ=$freq_normal;TUMOR_FREQ=$freq_tumor";
	}
	elsif($sFormat=~/DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50/) 
	{ 
		$varType="INDEL"; 
		$freq_normal=getINDELFreq($v_normal);
		$freq_tumor =getINDELFreq($v_tumor);
		$info="$info;NORMAL_FREQ=$freq_normal;TUMOR_FREQ=$freq_tumor";
	}
    }
    $F[7]=$info;
    print join("\t",@F)."\n";
    # For SNV:
    ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
    # ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
    # ##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
    # ##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
    # ##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
    # ##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
    # ##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
    # ##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
    # ##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
    #
    #1       1290016 .       T       A       .       PASS    NT=ref;QSS=74;QSS_NT=74;SGT=TT->AT;SOMATIC;TQSS=2;TQSS_NT=2     DP:FDP:SDP:SUBDP:AU:CU:GU:TU    129:0:0:0:0,0:0,0:1,1:128,128   99:0:0:0:18,18:1,1:1,1:79,81
    #1       1387689 .       G       A       .       PASS    NT=ref;QSS=73;QSS_NT=73;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1     DP:FDP:SDP:SUBDP:AU:CU:GU:TU    155:0:0:0:0,0:0,0:155,156:0,0   120:0:0:0:9,9:0,0:111,115:0,0
    #
    # For INDEL:
    ###FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1">
    ###FORMAT=<ID=DP2,Number=1,Type=Integer,Description="Read depth for tier2">
    ###FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2">
    ###FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">
    ###FORMAT=<ID=TOR,Number=2,Type=Integer,Description="Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2">
    ###FORMAT=<ID=DP50,Number=1,Type=Float,Description="Average tier1 read depth within 50 bases">
    ###FORMAT=<ID=FDP50,Number=1,Type=Float,Description="Average tier1 number of basecalls filtered from original read depth within 50 bases">
    ###FORMAT=<ID=SUBDP50,Number=1,Type=Float,Description="Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases">
    ###FILTER=<ID=Repeat,Description="Sequence repeat of more than 8x in the reference sequence">
    ###FILTER=<ID=iHpol,Description="Indel overlaps an interrupted homopolymer longer than 14x in the reference sequence">
    ###FILTER=<ID=BCNoise,Description="Average fraction of filtered basecalls within 50 bases of the indel exceeds 0.3">
    ###FILTER=<ID=QSI_ref,Description="Normal sample is not homozygous ref or sindel Q-score < 30, ie calls with NT!=ref or QSI_NT < 30">
    ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
    #1       120491661       .       ACT     A       .       PASS    IC=1;IHP=3;NT=ref;QSI=69;QSI_NT=69;RC=2;RU=CT;SGT=ref->het;SOMATIC;TQSI=1;TQSI_NT=1     DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50   155:155:156,158:0,0:3,3:151.35:0.00:0.00        212:212:174,174:23,23:13,13:202.02:0.00:0.00
    #2       179457265       .       C       CT      .       PASS    IC=7;IHP=8;NT=ref;QSI=308;QSI_NT=252;RC=6;RU=T;SGT=ref->het;SOMATIC;TQSI=2;TQSI_NT=1    DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50   424:424:397,398:0,0:22,22:428.93:0.43:0.00      454:454:297,299:93,93:63,63:454.42:0.00:0.00
    #
}
###################################################################

sub getSNVFreq
{
	my ($ref,$alt,$str)=@_;
	$alt=(split /,/,$alt)[0];
	my @k=split /:/,$str;
	my $refDepth=0;
	my $altDepth=0;
	if(   $ref eq "A") { $refDepth=(split /,/,$k[4])[1]; }
	elsif($ref eq "C") { $refDepth=(split /,/,$k[5])[1]; }
	elsif($ref eq "G") { $refDepth=(split /,/,$k[6])[1]; }
	elsif($ref eq "T") { $refDepth=(split /,/,$k[7])[1]; }

	if(   $alt eq "A") { $altDepth=(split /,/,$k[4])[1]; }
	elsif($alt eq "C") { $altDepth=(split /,/,$k[5])[1]; }
	elsif($alt eq "G") { $altDepth=(split /,/,$k[6])[1]; }
	elsif($alt eq "T") { $altDepth=(split /,/,$k[7])[1]; }

	return $altDepth/($refDepth+$altDepth);
}
sub getINDELFreq
{
	my ($str)=@_;
	my @k=split /:/,$str;
	my $refDepth=0;
	my $altDepth=0;
	$altDepth=(split /,/,$k[3])[1];
	$refDepth=$k[1]-$altDepth;
        #alt_dp=(l[10].split(':')[3]).split(',')[1]
        #ref_dp=str(int(l[10].split(':')[0])-int(alt_dp))

	return $altDepth/($refDepth+$altDepth);
}

sub usage
{
    die `pod2text $0`;
}
